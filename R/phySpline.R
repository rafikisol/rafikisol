#'@title analytically Solvable Physics-informed Categorical Depth Function
#'
#'@name phySpline
#'
#'@author Trevan Flynn
#'@description
#'Physics-informed spline that is fully analytically solvable and continuous in its 3 loss functions (data fidelity/mass-preservation, smoothing term, and physics term)
#'will interpolate "unseen" transition classes and not just the modal class. Informs class behavior and only assigns classes after all calculations have been performed.
#'Can run splines on a massive amount of soil data relatively fast and extract information that is scientific, not just algorithmic. !!BETA VERSION!!
#'
#'@param df a data frame with at least 4 columns with the first 3 named c("id","top","bottom"). The id column should have at least
#'@param class.var the character of the class to harmonise
#'@param lam the smoothing term [0, 1]
#'@param alpha the physics-informed term [0, 1]
#'@param d the vector to harmonise the soil depths to
#'@param order the latent ordering to map the classes to (default = class level order)
#'@param gthreshold the threshold of which the physics-informed loss function will start being applied
#'@param mode either "classification" or "continuous" for classification class.var must be a factor for now
#'
#'@return A list with 3 harmonised data frames with the classes, latent classes and uncertainties for each interval, a data frame of observations vs interpolations
#'and the lookup table of the data. It also returns the coefficients to be use in phySpline_Interp() as the function does not work off of 1 cm increments.
#'
#'@references Flynn, T., Rasaei, Z., Kostecki, R., (2025). Bayes’ Inference and Nearest Neighbor Splines to Harmonise Soil Texture Classes to International Depth Standards. https://doi.org/10.2139/ssrn.5078983
#'
#'@export
phySpline <- function(df, class.var, lam = 0.1, alpha = 0.1,
                      d = c(0,5,15,30,60,100,200),
                      order = levels(df[[class.var]]),
                      gthreshold = 30,
                      vlow = NULL, vhigh = NULL,
                      mode = c("classification", "continuous"),
                      show.progress = TRUE) {

  mode <- match.arg(mode)

  # --- 1. SETUP MAPPING & CONSTRAINTS ---
  if(mode == "classification"){
    if(!is.ordered(df[[class.var]])){
      df[[class.var]] <- factor(df[[class.var]], levels = order, ordered = TRUE)
    }
    df$post <- as.numeric(df[[class.var]])
    lookup <- data.frame(class = levels(df[[class.var]]),
                         numeric = seq_along(levels(df[[class.var]])))
    ycol <- "post"
  } else {
    ycol <- class.var
  }

  vlow_val <- if(!is.null(vlow)) vlow else min(df[[ycol]], na.rm = TRUE)
  vhigh_val <- if(!is.null(vhigh)) vhigh else max(df[[ycol]], na.rm = TRUE)

  profiles <- split(df, df[,1])
  n_profiles <- length(profiles)

  # Initialize storage
  coeff_list <- vector("list", n_profiles)
  obs_pred_list <- vector("list", n_profiles)

  if(show.progress) pb <- txtProgressBar(min=0, max=n_profiles, style=3)

  # --- 2. SOLVER LOOP ---
  for(p in seq_along(profiles)){
    prof <- profiles[[p]]
    u <- prof$top; v <- prof$bottom; y <- prof[[ycol]]
    n <- length(y)

    # Minimal data fallback
    if(n < 2){
      beta <- c(y[1], 0, 0)
    } else {
      # --- Analytical Fidelity Matrix ---
      A <- matrix(0, nrow = n, ncol = 3 * n)
      for(j in seq_len(n)){
        z1 <- u[j]; z2 <- v[j]; idx <- (3*(j-1)+1):(3*j)
        h_eff <- z2 - z1
        A[j, idx] <- c(1, (z2^2 - z1^2)/(2*h_eff), (z2^3 - z1^3)/(3*h_eff))
      }

      # Smoothness penalty
      R_mat <- matrix(0, 3*n, 3*n)
      for(j in seq_len(n)){
        h <- v[j]-u[j]; idx <- (3*j-2):(3*j)
        R_mat[idx, idx] <- matrix(c(0,0,0, 0,h,h^2, 0,h^2,4/3*h^3), nrow=3, byrow=TRUE)
      }

      # Continuity & physics
      G <- matrix(0, nrow = 2*(n-1), ncol = 3*n)
      for(j in 1:(n-1)){
        h_j <- v[j] - u[j]; idx_c <- (3*j-2):(3*j); idx_n <- (3*j+1):(3*j+3)
        w <- if(v[j]>=gthreshold) 1 else 0.1
        G[2*j-1, idx_c] <- c(1,h_j,h_j^2)*w; G[2*j-1, idx_n] <- c(-1,0,0)*w
        G[2*j, idx_c] <- c(0,1,2*h_j)*w; G[2*j, idx_n] <- c(0,-1,0)*w
      }

      scale_f <- sqrt(sum((t(A)%*%A)^2)) / max(1e-6, sqrt(sum((t(G)%*%G)^2)))
      G_s <- G*scale_f

      lhs <- t(A)%*%A + lam*R_mat + alpha*t(G_s)%*%G_s
      rhs <- t(A)%*%y

      # Gravitational / latent directional pull
      if(mode=="classification"){
        dy <- diff(y)
        L_pull <- rep(0, 3*n)
        for(j in 1:(n-1)){
          L_pull[3*j-1] <- dy[j]; L_pull[3*j] <- dy[j]/2
        }
        rhs <- rhs + alpha*L_pull
      }

      beta <- tryCatch(solve(lhs, rhs), error=function(e) rep(NA,3*n))
    }

    # Store coefficients
    coeff_list[[p]] <- data.frame(alpha=beta[seq(1,3*n,3)],
                                  b=beta[seq(2,3*n,3)],
                                  gamma=beta[seq(3,3*n,3)],
                                  top=u, bottom=v)

    # Observed vs predicted (latent or class)
    pred_latent <- as.numeric(A %*% beta[1:(3*n)])
    if(mode=="classification"){
      pred_class <- sapply(pred_latent, function(val){
        if(is.na(val)) return(NA)
        lookup$class[which.min(abs(val - lookup$numeric))]
      })
      obs_pred_list[[p]] <- data.frame(id=prof[[1]][1],
                                       top=u, bottom=v,
                                       obs=prof[[class.var]],
                                       pred=pred_class)
    } else {
      obs_pred_list[[p]] <- data.frame(id=prof[[1]][1],
                                       top=u, bottom=v,
                                       obs=y, pred=pred_latent)
    }

    if(show.progress) setTxtProgressBar(pb,p)
  }
  if(show.progress) close(pb)

  # --- 3. HARMONIZATION VIA DEFINITE INTEGRALS ---
  harmonized_df <- data.frame(id=names(profiles))
  harmonized_raw <- data.frame(id=names(profiles))
  harmonized_uncert_df <- data.frame(id=names(profiles))

  for(i in 1:(length(d)-1)){
    z1 <- d[i]; z2 <- d[i+1]; col_name <- paste0(z1,"-",z2,"cm")

    harmonized_raw[[col_name]] <- sapply(seq_along(coeff_list), function(p){
      cf <- coeff_list[[p]]; p_max <- max(cf$bottom)
      if(z1>=p_max) return(NA)
      z2_eff <- min(z2,p_max); total_area <- 0
      for(j in 1:nrow(cf)){
        h_s <- max(cf$top[j],z1); h_e <- min(cf$bottom[j],z2_eff)
        if(h_e>h_s){
          x1 <- h_s - cf$top[j]; x2 <- h_e - cf$top[j]
          area <- cf$alpha[j]*(x2-x1) + (cf$b[j]/2)*(x2^2-x1^2) + (cf$gamma[j]/3)*(x2^3-x1^3)
          total_area <- total_area + area
        }
      }
      total_area / (z2_eff-z1)
    })

    if(mode=="classification"){
      harmonized_df[[col_name]] <- sapply(harmonized_raw[[col_name]], function(val){
        if(is.na(val)) return(NA)
        lookup$class[which.min(abs(val - lookup$numeric))]
      })
      harmonized_uncert_df[[col_name]] <- sapply(harmonized_raw[[col_name]], function(val){
        if(is.na(val)) return(NA); min(abs(val - lookup$numeric))
      })
    } else {
      harmonized_df[[col_name]] <- harmonized_raw[[col_name]]
    }
  }

  return(list(
    harmonised = harmonized_df,
    harmonised_raw = harmonized_raw,
    uncertainty = harmonized_uncert_df,
    obs.preds = do.call(rbind, obs_pred_list),
    lookup = if(mode=="classification") lookup else NULL,
    coeffs = coeff_list
  ))
}
