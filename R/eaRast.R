#'@title A method to run equal-area splines on a stack or rasters
#'
#'@name eaRast
#'
#'@description
#'This algorithm takes a stack of rasters or spatrasters (decreasing depths),
#'and runs an ea_spline (Malone et al., 2019) at each pixel. To get the final
#'value it then integrates the spline (e.g., carbon stock), which is much more
#'accurate than summing 1 cm intervals.
#'
#'@param x Is the stack of spatrasters to be used to reduce.
#'@param spc A soilProfileCollection created from the spatraster to reduce.
#'@param depths Is a vector of the depths (cm) to integrate (between 1:200).
#'@param var Character of what property the spline should be run on.
#'@param type A character specifying the type of calculations to perform on the spline.
#'can either be "integrate" (cumulative sum), "sum" (descrete sum), "mean",
#'or "median".
#'
#'@return A list consisting of the reduced spatraster and all the elements
#'produced from the ea_spline function from the ithir R package (Malone et al., 2009).
#'
#'@references https://10.1016/S0016-7061(99)00003-8
#'@references http://dx.doi.org/10.1016/j.geoderma.2009.10.007.
#'
#'@export

eaRast = function(x, spc, depths, var, type = 'integrate'){

  #get coordinates
  coord = na.omit(as.data.frame(x, xy = T,)[, 1:2])

  #run spline
  spl = ithir::ea_spline(spc, var.name = var, vlow = 0)

  #get predictions
  spred = as.data.frame(spl$var.1cm)

  #depth matrix
  dp = matrix(NA, nrow = nrow(spred), ncol = ncol(spred))
  dp = as.data.frame(apply(dp, 2, function(x) 1:200)) #add values

  #integrate
  if(type == "integrate"){
    #integrate columns
    for(i in 1:ncol(spred)){

      coord[i, 3] = integrate(approxfun(spred[depths, i], dp[depths, i]), lower = min(spred[depths, i]), upper = max(spred[depths, i]),
                              subdivisions = max(dp[depths, i]), stop.on.error = F)$value
    }
    #make spatraster and adjust for depth
    ras = rast(coord, crs = crs(x))*max(depths)/100
    #resample to original
    res = resample(ras, x)
    #fill voids
    ras = focal(ras, w = 9, fun = "median", na.policy = "only", fillvalue = NA)

    return(setNames(list(ras, spl), c("predictions", "models")))
  }

  #descete sum
  if(type == "sum"){

    for(i in 1:ncol(spred)){
      coord[i, 3] = sum(spred[depths, i])
    }

    ras = rast(coord, crs = crs(x))*max(depths)/100
    res = resample(ras, x)
    ras = focal(ras, w = 9, fun = "median", na.policy = "only")

    return(setNames(list(ras, spl), c("predictions", "models")))
  }

  #mean
  if(type == "mean"){
    #loop around values
    for(i in 1:ncol(spred)){
      coord[i, 3] = mean(spred[depths, i])
    }
    #make spatraster and adjust to depth
    ras = rast(coord, crs = crs(x))*max(depths)/100
    #resample to original
    res = resample(ras, x)
    #fill voids
    ras = focal(ras, w = 9, fun = "median", na.policy = "only", fillvalue = NA)

    return(setNames(list(ras, spl), c("predictions", "models")))
  }

  #median
  if(type == "median"){

    #fill in correct order
    for(i in 1:ncol(spred)){
      coord[i, 3] = median(spred[depths, i])
    }

    #make into spatraster
    ras = rast(coord, crs = crs(x))*max(depths)/100
    res = resample(ras, x)

    #fill voids
    ras = focal(ras, w = 9, fun = "median", na.policy = "only", fillvalue = NA)
    return(setNames(list(ras, spl), c("predictions", "models")))
  }
}
#End
