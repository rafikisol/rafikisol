#'@title analytically Solvable Physics-informed Categorical Depth Function
#'
#'@name phySpline
#'
#'@author Trevan Flynn
#'@description
#'Physics-informed spline that is fully analytically solvable and infinitely differentiable. The spline self-regulates of the cumulative mass of the overlaying horizons
#'stiffening the spline with depth yet the lambda has a large influence and can never over constrain the spline. The spline fits the values by minimising the jacobian
#'energy held in the soil vertical gradients giving meaning to the curve. Since everything is a continuous function its extremely efficient and difficult to over constrain
#'however, it is still in the beta version.
#'
#'@param df a data frame with at least 4 columns with the first 3 named c("id","top","bottom"). The id column should have at least
#'@param class.var the character of the class to harmonise
#'@param lam the smoothing term [0, 1]
#'@param d the vector of soil depths intervals to harmonise to.
#'@param order the latent ordering to map the classes to (default = class level order)
#'@param mode either "classification" or "continuous" for classification class.var must be a factor for now
#'
#'@return A list with 3 harmonised data frames with the classes, latent classes and uncertainties for each interval, a data frame of observations vs interpolations
#'and the lookup table of the data. It also returns the coefficients to be use in phySpline_Interp() as the function does not work off of 1 cm increments.
#'
#'@references Flynn, T., Rasaei, Z., Kostecki, R., (2025). Bayes’ Inference and Nearest Neighbor Splines to Harmonise Soil Texture Classes to International Depth Standards. https://doi.org/10.2139/ssrn.5078983
#'
#'@export
phySpline <- function(df,
                           class.var,
                           lam = 0.1,
                           d = c(0,5,15,30,60,100,200),
                           mode = c("classification", "continuous"),
                           order = levels(df[[class.var]]),
                           show.progress = TRUE) {

  mode <- match.arg(mode)

  # ------------------------------------------------------------
  # 1. Classification setup
  # ------------------------------------------------------------
  if (mode == "classification") {
    if (!is.ordered(df[[class.var]])) {
      df[[class.var]] <- factor(df[[class.var]], levels = order, ordered = TRUE)
    }
    df$post <- as.numeric(df[[class.var]])
    ycol <- "post"
    lookup <- data.frame(
      class   = levels(df[[class.var]]),
      numeric = seq_along(levels(df[[class.var]]))
    )
  } else {
    ycol <- class.var
    lookup <- NULL
  }

  profiles <- split(df, df[,1])
  n_profiles <- length(profiles)

  # ------------------------------------------------------------
  # 2. Storage
  # ------------------------------------------------------------
  coeff_list      <- vector("list", n_profiles)
  full_coeff_list <- vector("list", n_profiles)
  obs_pred_list   <- vector("list", n_profiles)
  diag_list       <- vector("list", n_profiles)

  fit_errors <- data.frame(
    id = names(profiles),
    rmse = NA_real_,
    rmse_norm = NA_real_
  )

  if (show.progress) pb <- txtProgressBar(0, n_profiles, style = 3)

  # ------------------------------------------------------------
  # 3. Profile loop
  # ------------------------------------------------------------
  for (p in seq_along(profiles)) {

    prof <- profiles[[p]][order(profiles[[p]]$top), ]
    prof$is_ghost <- FALSE

    # ---- Insert ghosts
    filled <- list()
    if (prof$top[1] > 0) {
      g <- prof[1, ]; g$top <- 0; g$bottom <- prof$top[1]
      g[[ycol]] <- prof[[ycol]][1]; g$is_ghost <- TRUE
      filled[[1]] <- g
    }
    for (i in seq_len(nrow(prof))) {
      filled[[length(filled)+1]] <- prof[i, ]
      if (i < nrow(prof) && prof$bottom[i] < prof$top[i+1]) {
        g <- prof[i, ]; g$top <- prof$bottom[i]; g$bottom <- prof$top[i+1]
        g[[ycol]] <- median(prof[[ycol]][c(i,i+1)], na.rm=TRUE); g$is_ghost <- TRUE
        filled[[length(filled)+1]] <- g
      }
    }
    prof <- do.call(rbind, filled)

    # ---- Geometry
    u  <- prof$top; v <- prof$bottom; y <- prof[[ycol]]
    hi <- v - u; H <- sum(hi); n <- length(y)

    # ---- Primary: analytical solver with fallback
    if (n > 1) {
      # 1. Geometry Matrices (Already defined u, v, hi, H, n above)
      A <- matrix(0, n, 3*n)
      for (j in seq_len(n)) {
        h <- hi[j]
        A[j, (3*j-2):(3*j)] <- c(1, h/2, h^2/3)
      }

      # Function for Solver: Define the Refined Euler Resistance Process (increases speed)
      calc_stiffness <- function(z_top, z_bot, k) {
        h_norm <- z_bot - z_top
        # Numerical safety: If horizon is ultra-thin, return the point-stiffness
        if (h_norm < 1e-6) return(exp(k * z_top))
        # The Integral Mean: Captures the total resistance across the layer
        return((exp(k * z_bot) - exp(k * z_top)) / (k * h_norm))
      }

      # 2. Inside the Profile Loop (after geometry setup)
      z_bot_norm <- cumsum(hi) / H
      z_top_norm <- c(0, z_bot_norm[-n])
      k <- exp(1)
      
      # regularisation matrix - collapses to Jacobian energy holomorphic integral
      R <- matrix(0, 3*n, 3*n)
      
      for(j in seq_len(n)) {
        idx <- (3*j-2):(3*j)
        h <- hi[j]
        # (Penalizes first derivative weighted by squared stiffness)
        R_block <- matrix(c(0, 0, 0,
                            0, h, h^2,
                            0, h^2, 4/3*h^3), 3, 3, byrow = TRUE)
        # calculate the integral on the integral collapsed on the inginite decay
        R[idx, idx] <- calc_stiffness(z_top_norm[j], z_bot_norm[j], k)^2 * R_block
      }
      # Continuity constraints
      Gc <- matrix(0, 2*(n-1), 3*n)
      for (j in 1:(n-1)) {
        h <- hi[j]; i1 <- (3*j-2):(3*j); i2 <- (3*j+1):(3*j+3)
        Gc[2*j-1, i1] <- c(1,h,h^2); Gc[2*j-1, i2] <- c(-1,0,0)
        Gc[2*j,   i1] <- c(0,1,2*h); Gc[2*j, i2] <- c(0,-1,0)
      }
      Gtop <- matrix(0,1,3*n); Gtop[1,2] <- 1
      Gbot <- matrix(0,1,3*n); Gbot[1,(3*n-1):(3*n)] <- c(1,2*hi[n])
      G <- rbind(Gc,Gtop,Gbot)

      # Scaling
      sA <- sqrt(sum((t(A)%*%A)^2))
      R_s <- R*(sA/max(1e-8,sqrt(sum(R^2))))
      G_s <- G*(sA/max(1e-8,sqrt(sum((t(G)%*%G)^2))))

      # Null-space
      qrG <- qr(t(G_s), tol=1e-7)
      Q <- qr.Q(qrG, complete=TRUE)
      N <- Q[, (qrG$rank+1):ncol(Q), drop=FALSE]

      # ---- Solve with fallback
      solve_beta <- function() {
        if (ncol(N) > 0) {
          AN <- A %*% N
          lhs <- t(AN)%*%AN/n + lam*t(N)%*%R_s%*%N
          rhs <- t(AN)%*%y/n
          N %*% tryCatch(solve(lhs,rhs), error=function(e) qr.coef(qr(lhs),rhs))
        } else rep(0,3*n)
      }

      beta <- tryCatch({
        solve_beta()
      }, error=function(e){
        # fallback: simple numeric integration
        b <- numeric(3*n)
        for(j in 1:n) b[(3*j-2):(3*j)] <- c(y[j],0,0)
        b
      })
      # tells if we couldnt solve analytically                   
      method_used <- if(exists("solve_beta")) "analytical" else "numeric"

      n_segments <- length(beta) / 3
      # buld coefficient matrix                   
      cf_all <- data.frame(
        alpha  = beta[seq(1, 3 * n_segments, 3)],
        b      = beta[seq(2, 3 * n_segments, 3)],
        gamma  = beta[seq(3, 3 * n_segments, 3)],
        top    = u[1:n_segments],    # Crucial: slice to match n_segments
        bottom = v[1:n_segments]     # Crucial: slice to match n_segments
      )
      # fit the predictions from coefficients
      pred <- cf_all$alpha + cf_all$b*hi/2 + cf_all$gamma*hi^2/3

    } else {
      # if just one horizon
      cf_all <- data.frame(alpha=y,b=0,gamma=0,top=u,bottom=v,is_ghost=FALSE)
      pred <- y; beta <- rep(0,3*n)
      method_used <- "single"
    }

    # ---- Store
    full_coeff_list[[p]] <- cf_all
    coeff_list[[p]] <- cf_all

    # Get RMSE of residuals                    
    real <- !prof$is_ghost
    rmse <- sqrt(mean((y[real]-pred[real])^2, na.rm=TRUE))
    fit_errors$rmse[p] <- rmse
    fit_errors$rmse_norm[p] <- rmse / max(diff(range(y[real])), 1e-10)

    # Classification
    if (mode=="classification") {
      pred_class <- lookup$class[
        max.col(-abs(outer(pred, lookup$numeric,"-")), ties.method="first")
      ]
    }
    # fitted values with latent residuals                     
    obs_pred_list[[p]] <- data.frame(
      id = names(profiles)[p],
      top = u[real],
      bottom = v[real],
      y = y[real], # want exact latent values
      obs = if(mode=="classification") prof[[class.var]][real] else y[real],
      pred = if(mode=="classification") pred_class[real] else pred[real],
      value = pred[real],
      resid = y[real]-pred[real]
    )
    # Tells if analytical, over-constrained and errors
    diag_list[[p]] <- data.frame(
      id = names(profiles)[p],
      method_used = method_used,
      max_constraint = max(abs(G_s %*% beta)),
      rmse = rmse
    )

    if(show.progress) setTxtProgressBar(pb,p)
  }

  if(show.progress) close(pb)

  # ------------------------------------------------------------
  # 4. Harmonisation - integrate to get values of of classes and uncertainty
  # ------------------------------------------------------------
  harmonised <- data.frame(id=names(profiles))
  uncert <- data.frame(id=names(profiles))
  
  for(i in seq_len(length(d)-1)){
    z1 <- d[i]; z2 <- d[i+1]; nm <- paste0(z1,"-",z2,"cm")

    vals <- sapply(full_coeff_list, function(cf){
      if(any(is.na(cf$alpha))) return(NA)
      area <- 0
      for(j in seq_len(nrow(cf))){
        a <- max(cf$top[j], z1); b <- min(cf$bottom[j], z2)
        if(b>a){
          x1 <- a - cf$top[j]; x2 <- b - cf$top[j]
          area <- area + cf$alpha[j]*(x2-x1) + cf$b[j]/2*(x2^2-x1^2) + cf$gamma[j]/3*(x2^3-x1^3)
        }
      }
      area/(z2-z1)
    })

    if(mode=="classification"){
      closest <- lookup$numeric[max.col(-abs(outer(vals, lookup$numeric,"-")))]
      harmonised[[nm]] <- factor(lookup$class[match(closest, lookup$numeric)], levels=lookup$class)
      uncert[[nm]] <- pmin(1, abs(vals-closest))
    } else {
      harmonised[[nm]] <- vals
    }
  }
  # add soil depth to harmonised data frames
  harmonised$soil_depth <- sapply(coeff_list, function(cf) max(cf$bottom))
  uncert$soil_depth <- harmonised$soil_depth

  # ------------------------------------------------------------
  # 5. Return
  # ------------------------------------------------------------
  list(
    harmonised = harmonised,
    uncertainty = uncert,
    obs.preds = do.call(rbind, obs_pred_list),
    lookup = lookup,
    coeffs = coeff_list,
    errors = fit_errors,
    diagnostics = do.call(rbind, diag_list)
  )
}
