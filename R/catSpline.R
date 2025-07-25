#'@title Function for categorical splines
#'@name catSpline
#'@author Trevan Flynn
#'
#'@description Function runs equal-area splines on the probabilities of soil categorical data
#'and then classifies the soil class through a 1D nearest neighbor approach. Based of
#'the code written by Brendon Melone, catSpline is only for categorical data and
#'uses the softmedian for the harmonised depths to maintain the differential qualities
#'while aligning with the modal for categorical data.
#'
#'@param obj is a data frame with columns of profile top of horizon and bottom of horizon
#'These columns must be named id, top and bottom respectively to be turned into a spc.f
#'@param class.var is the name of the factor to run the spline.
#'@param type type of output being "class" with the classification or "probs" with the raw probabilities.
#'@param lam is the smoothing factor of the spline
#'@param beta is how smooth the softmedian value is where lower is smoother and
#'higher is closer to the true median.
#'@param d is the depths to harmonise to.
#'@param vlow is the lower boundary value (lowest = 0).
#'@param vhigh is the higher boundary value (highest = 0).
#'
#'@returns a list of harmonised posterior probabilites, observed vs predictions,
#'predictions in 1 cm increments (for figures) and the lookup table with probabilities.
#'
#'@references Malone, B.P., McBratney, A.B., Minasny, B., Laslett, G.M. (2009) Mapping continuous depth functions of soil carbon storage and available water capacity. Geoderma, 154(1-2): 138-152.
#'@references Bishop, T.F.A., McBratney, A.B., Laslett, G.M., (1999) Modelling soil attribute depth functions with equal-area quadratic smoothing splines. Geoderma, 91(1-2): 27-45.
#'
#'@export
# Spline fitting for horizon data (Matlab Code converted to R by Brendan Malone)
catSpline<- function(obj, class.var, type = "class", beta = 1, lam = 0.1, d = c(0,5,15,30,60,100,200), vlow = 0, vhigh = 1,show.progress=TRUE){

  #load aqp -- helps keep it consistant
  require(aqp)

  #calculate differential median
  softmedian <- function(x, beta = 1) {
    # x: Numeric vector
    # beta: Softness parameter (higher beta = closer to true median)
    weights <- exp(-beta * abs(x - median(x)))
    sum(weights * x) / sum(weights)
  }

  #1D nearest neibor classification
  closest_value <- function(x, values) {
    if(is.numeric(x)){
      values[which.min(abs(values - x))]
    }
  }

  second_value <- function(x, values) {
    # Calculate absolute differences
    diffs <- abs(values - x)

    # Find the indices of the two smallest differences
    closest_indices <- order(diffs)[1:2]

    # Return the second closest value
    x[closest_indices[2]]
  }

  ###############################################################################
  #get prior probs
  tab = table(obj[[class.var]])
  prob = tab/sum(tab)
  loglik = tab * log(prob)
  post = prob * loglik
  post = post/sum(post)
  lookup = as.data.frame(cbind(names(tab), post))
  lookup$post = as.numeric(lookup$post)
  names(lookup) = c(class.var, "post")

  var.name = "post"

  obj = merge(obj, lookup, by = class.var)

  aqp::depths(obj) = id ~ top + bottom
  ##############################################################################
  if (is(obj,"SoilProfileCollection") == TRUE){
    depthcols = obj@depthcols
    idcol = obj@idcol
    obj@horizons = obj@horizons[,c(idcol, depthcols, var.name)]
    d.dat<- as.data.frame(obj@horizons)}

  if (is(obj,"data.frame") == FALSE & is(obj,"SoilProfileCollection") == FALSE){
    stop("ERROR:: Data must be either class data.frame or SoilProfileCollection")}

  mxd<- max(d)
  sp_dat<-split(d.dat,d.dat[,1])

  # matrix of the continous splines for each data point
  m_fyfit<- matrix(NA,ncol=length(c(1:mxd)),nrow=length(sp_dat))

  # Matrix in which the averaged values of the spline are fitted. The depths are specified in the (d) object
  yave<- matrix(NA,ncol=length(d),nrow=length(sp_dat))

  # Matrix in which the sum of square errors of each lamda iteration for the working profile are stored
  sse<- matrix(NA,ncol=length(lam),nrow=1)

  # Matrix in which the sum of square errors for eac h lambda iteration for each profile are stored
  sset<- matrix(NA,ncol=2,nrow=length(sp_dat))

  #Profile ids
  mat_id<- d.dat[0,]

  #combined data frame for observations and spline predictions
  dave<- d.dat[1,]
  dave$predicted<- 0
  dave$FID<- 0
  dave<- dave[0,]

  ## Fit splines profile by profile:
  if (show.progress) pb <- utils::txtProgressBar(min=0, max=length(sp_dat), style=3)
  cnt<- 1

  for(st in 1:length(sp_dat)){
    subs<-sp_dat[[st]]  # subset the profile required
    subs<-as.matrix(subs)
    mat_id[st,1]<- subs[1,1]


    # manipulate the profile data to the required form
    ir<- c(1:nrow(subs))
    ir<-as.matrix(t(ir))
    u<- as.numeric(subs[ir,2])
    u<-as.matrix(t(u))   # upper
    v<- as.numeric(subs[ir,3])
    v<-as.matrix(t(v))   # lower
    y<- as.numeric(subs[ir,4])
    y<-as.matrix(t(y))   # concentration
    n<- length(y);       # number of observations in the profile
    d<- t(d)


    ############################################################################################################################################################
    ### routine for handling profiles with one observation
    if (n == 1)
    {dave[cnt:(cnt-1+nrow(subs)),1:4]<- subs
    dave[cnt:(cnt-1+nrow(subs)),5]<- y
    dave[cnt:(cnt-1+nrow(subs)),6]<- st
    xfit<- as.matrix(t(c(1:mxd))) #spline will be interpolated onto these depths (1cm res)
    nj<- max(v)
    if (nj > mxd)
    {nj<- mxd}
    yfit<- xfit
    yfit[,1:nj]<- y   #values extrapolated onto yfit
    if (nj < mxd)
    {yfit[,(nj+1):mxd]=-9999}
    m_fyfit[st,]<- yfit


    # Averages of the spline at specified depths
    nd<- length(d)-1  # number of depth intervals
    dl<-d+1     #  increase d by 1

    for(cj in 1:nd){
      xd1<- dl[cj]
      xd2<- dl[cj+1]-1
      if (nj>=xd1 & nj<=xd2)
      {xd2<- nj-1
      yave[st,cj]<- softmedian(yfit[,xd1:xd2])}
      else
      {yave[st,cj]<- softmedian(yfit[,xd1:xd2])}   #median of the spline at the specified depth intervals
      yave[st,cj+1]<- max(v)#maximum soil depth
    }

    cnt<- cnt+nrow(subs)

    # error
    sset[st,1:2] <- 0
    }

    ###############################################################################################################################################################

    ##Profile set up
    else{
      dave[cnt:(cnt-1+nrow(subs)),1:4]<- subs
      dave[cnt:(cnt-1+nrow(subs)),6]<- st
      ## ESTIMATION OF SPLINE PARAMETERS
      np1 <- n+1  # number of interval boundaries
      nm1 <- n-1
      delta <- v-u  # depths of each layer
      del <- c(u[2:n],u[n])-v   # del is (u1-v0,u2-v1, ...)

      ## create the (n-1)x(n-1) matrix r; first create r with 1's on the diagonal and upper diagonal, and 0's elsewhere
      r <- matrix(0,ncol=nm1,nrow=nm1)

      for(dig in 1:nm1){
        r[dig,dig]<-1
      }

      for(udig in 1:nm1-1){
        r[udig,udig+1]<-1
      }

      ## then create a diagonal matrix d2 of differences to premultiply the current r
      d2 <- matrix(0, ncol=nm1, nrow=nm1)
      diag(d2) <- delta[2:n]  # delta = depth of each layer

      ## then premultiply and add the transpose; this gives half of r
      r <- d2 %*% r
      r <- r + t(r)

      ## then create a new diagonal matrix for differences to add to the diagonal
      d1 <- matrix(0, ncol=nm1, nrow=nm1)
      diag(d1) <- delta[1:nm1]  # delta = depth of each layer

      d3 <- matrix(0, ncol=nm1, nrow=nm1)
      diag(d3) <- del[1:nm1]  # del =  differences

      r <- r+2*d1 + 6*d3

      ## create the (n-1)xn matrix q
      q <- matrix(0,ncol=n,nrow=n)

      for(dig in 1:n){
        q[dig,dig]<- -1
      }

      for(udig in 1:n-1){
        q[udig,udig+1]<-1
      }

      q <- q[1:nm1,1:n]
      dim.mat <- matrix(q[],ncol=length(1:n),nrow=length(1:nm1))

      ## inverse of r
      rinv <- try(solve(r), TRUE)

      ## Note: in same cases this will fail due to singular matrix problems, hence you need to check if the object is meaningfull:
      if(is.matrix(rinv)){
        ## identity matrix i
        ind <- diag(n)

        #create the matrix coefficent z
        pr.mat <- matrix(0,ncol=length(1:nm1),nrow=length(1:n))
        pr.mat[] <- 6*n*lam
        fdub <- pr.mat*t(dim.mat)%*%rinv
        z <- fdub%*%dim.mat+ind

        #solve for the fitted layer means
        sbar <- solve(z,t(y))
        dave[cnt:(cnt-1+nrow(subs)),5]<- sbar
        cnt<- cnt+nrow(subs)


        #calculate the fitted value at the knots
        b <- 6*rinv%*%dim.mat%*% sbar
        b0 <- rbind(0,b) # add a row to top = 0
        b1 <- rbind(b,0) # add a row to bottom = 0
        gamma <- (b1-b0) / t(2*delta)
        alfa <- sbar-b0 * t(delta) / 2-gamma * t(delta)^2/3

        ###############################################################################################################################################################
        #Spline fit
        xfit<- as.matrix(t(c(1:mxd))) #1 cm estimates
        nj<- max(v)
        if (nj > mxd)
        {nj<- mxd}
        yfit<- xfit

        for(k in 1:nj){
          xd<-xfit[k]
          if (xd< u[1])
          {p<- alfa[1]} else
          {for (its in 1:n) {
            if(its < n)
            {tf2=as.numeric(xd>v[its] & xd<u[its+1])} else {tf2<-0}
            if (xd>=u[its] & xd<=v[its])
            {p=alfa[its]+b0[its]*(xd-u[its])+gamma[its]*(xd-u[its])^2} else if (tf2)
            {phi=alfa[its+1]-b1[its]*(u[its+1]-v[its])
            p=phi+b1[its]*(xd-v[its])}
          }}
          yfit[k]=p
        }

        if(nj < mxd)
        {yfit[,(nj+1):mxd]=NA}

        yfit[which(yfit > vhigh)] <- vhigh
        yfit[which(yfit < vlow)]  <-vlow
        m_fyfit[st,]<- yfit

        ## Averages of the spline at specified depths
        nd<- length(d)-1  # number of depth intervals
        dl<-d+1     #  increase d by 1

        for(cj in 1:nd){
          xd1<- dl[cj]
          xd2<- dl[cj+1]-1

          if(nj>=xd1 & nj<=xd2){
            xd2<- nj-1
            yave[st,cj]<- softmedian(yfit[,xd1:xd2])
          }
          else{
            yave[st,cj]<- softmedian(yfit[,xd1:xd2])} #softmedian of the spline at the specified depth intervals
          yave[st,cj+1]<- max(v)#maximum soil depth
        }

        ##Keep the original error of probabilities
        rmse <- sqrt(sum((t(y)-sbar)^2)/n)
        rmseiqr <- rmse/stats::IQR(y)
        sset[st,1] <- rmse
        sset[st,2] <- rmseiqr
      }
    }

    if(show.progress){utils::setTxtProgressBar(pb, st)}
  }
  if(show.progress){
    close(pb)
  }

  ## asthetics for output-----------------------------------------------------
  ## yave
  yave<- as.data.frame(yave)
  jmat<- matrix(NA,ncol=1,nrow=length(d))
  for (i in 1:length(d)-1) {
    a1<-paste(d[i],d[i+1],sep="-")
    a1<-paste(a1,"cm",sep=" ")
    jmat[i]<- a1}

  for (jj in 1:length(jmat)){
    names(yave)[jj]<- jmat[jj]
  }
  jmat[length(d)]<- "soil_depth"
  for (jj in 1:length(jmat)){
    names(yave)[jj]<- jmat[jj]
  }
  yave<- cbind(mat_id[,1],yave)
  names(yave)[1]<- "id"

  # Convert 'id' and 'soil_depth' to factors
  yave[, "id"] <- as.factor(yave[, "id"])
  yave[, "soil_depth"] <- as.factor(yave[, "soil_depth"])

  numColNames <- names(yave)[sapply(yave, is.numeric)]
  uncert <- yave

  # Apply closest_value to each numerical column and get uncertainties
  for (col_name in numColNames) {
    yave[[col_name]] <- sapply(yave[[col_name]], closest_value,
                               values = lookup$post,
                               USE.NAMES = FALSE)
  }

  for (col_name in numColNames) {
    # mapply returns vector or list; wrap in unlist() to flatten
    result <- mapply(function(x, y) {
      val <- tryCatch(abs(x - y), error = function(e) NA)
      if (length(val) == 0 || is.null(val)) NA else val
    }, x = uncert[[col_name]], y = yave[[col_name]])

    uncert[[col_name]] <- result
  }
  #if probabilities desired-------------------------------------------------------
  if(type == "probs"){
    # sset
    sset<- as.data.frame(sset)
    names(sset)<- c("rmse", "rmseiqr")

    # save outputs
    retval <- list(harmonised=yave,
                   uncertainty = uncert,
                   obs.preds=dave,
                   splineFitError=sset ,
                   var.1cm=t(m_fyfit))
    return(retval)}


  if(type == "class"){
    # Identify numerical columns by name

    # Convert the numerical columns to factors based on lookup$Class
    for (col_name in numColNames) {
      yave[[col_name]] <- factor(lookup$Class[match(yave[[col_name]], lookup$post)],
                                 levels = lookup$Class) # Explicitly set levels for consistent factors
    }
    #dave
    dave
    #Match probabilities to class
    closest_obvs <- sapply(dave$post, function(x) closest_value(as.numeric(x), lookup$post))
    dave$Obvs <- as.factor(lookup$Class[match(closest_obvs, lookup$post)])

    closest_pred <- sapply(dave$predicted, function(x) closest_value(x, lookup$post))
    dave$Pred <- factor(lookup$Class[match(closest_pred, lookup$post)],
                        levels = levels(dave$Obvs))

    dave = dave[, c("id", "top", "bottom", "post", "predicted", "Obvs", "Pred")]

    #1cm preds
    names(m_fyfit) = gsub("V", "", names(m_fyfit))

    closest_pred <- sapply(m_fyfit, function(x) closest_value(x, lookup$post))
    m_fyfit <- factor(lookup$Class[match(closest_pred, lookup$post)])

    #sset
    sset<- as.data.frame(sset)
    names(sset)<- c("rmse", "rmseiqr")

    # save outputs
    retval <- list(harmonised=yave,
                   uncertainty = uncert,
                   obs.preds = dave ,
                   var.1cm= t(m_fyfit),
                   lookup = lookup,
                   Error=sset)
    return(retval)
  }
}
