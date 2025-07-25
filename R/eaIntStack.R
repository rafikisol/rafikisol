#'@title A function to integrate equal area splines in a raster stack
#'
#'@name eaIntStack
#'
#'@author Trevan Flynn
#'
#'@description A function that takes a stack of rasters with the same properties
#'and reduces into an image through an equal area spline and takes the integral of
#'the results at each pixel.Therefore, its a depth function for a raster stack. Essentially it gives the total
#'for all depths combined in a precise manner and it is  more accurate than taking averages
#'or summation of the means.Created by Ponce-Hernandez et al. (1986), Bishop et al. (1999),  and implemented
#'and improved by Malone et al. (2009) where this function came from.
#'
#'@param stack Is a spatraster with each band of the same property with depth.
#'@param label A character of what is to be predicted.
#'@param depths A vector from the upper limit to lower limit of what should be integrated.
#'@param n Is the number of rasters in the stack.
#'@param upper Is a vector of the top depth of each measurement.
#'@param lower Is a vector of the lower depth of each measurement.
#'@param ... Is any additional arguments past to ea_spline() function
#'
#'@return Is an image with the integrated value from all images.
#'
#'@references  R. Ponce-Hernandez, F.H.C. Marriott, P.H.T. Beckett (1986).An improved method for reconstructing a soil profile from analyses of a small number of samples.
#'Journal of Soil Science, 37, pp. 455-467,https://doi.org/10.1111/j.1365-2389.1986.tb00377.x
#'@references T.F.A. Bishop, A.B. McBratney, G.M. Laslett, (1999). Modelling soil attribute depth functions with equal-area quadratic smoothing splines.
#'Geoderma, 91, pp. 27-45,https://10.1016/S0016-7061(99)00003-8
#'@references B.P. Malone, A.B. McBratney, B. Minasny, G.M. Laslett, (2009). Mapping continuous depth functions of soil carbon storage and available water capacity.
#'Geoderma, 154, pp. 138-152,https://10.1016/j.geoderma.2009.10.007
#'
#'@note Since it runs a spline and takes the integral for each pixel. This can take
#'some time for large or high resolution rasters. We recommend preforming the
#'function on a coarse image and downscaling through dissever() or rafikiReg().Or,
#'one could sample the raster and predict the total (e.g., carbon stock) over
#'covariates.
#'
#'@export


eaIntStack = function(stack, label, depths = 0:200, model = "ea", n,
                       upper = c(0, 5, 15, 30, 60, 100),
                       lower = c(5, 15, 30, 60, 100, 200), ...){

  #make an with unique values at each pixel for the profile variable
  r = stack[[1]] #make raster
  terra::values(r) = 1:na.omit(terra::values(stack[[1]])) #replace values with numbers
  names(r) = "profile" #name "profile" is a must

  #convert raster
  dat = as.data.frame(c(stack, r), xy = T)%>%
    dplyr::pivot_longer(., col = c(3:n), values_to = label)%>% #pivot longer
    na.omit()

  #put in depths - should equal the number of rasters in the stack
  dat$top = rep_len(upper, length.out = nrow(dat))
  dat$bottom = rep_len(lower, length.out = nrow(dat))

  if(model == "ea"){#change to profile collection

    aqp::depths(dat) = profile ~ top + bottom #spc
    aqp::site(dat) = ~ x + y #make sure each pixel is its own site
    aqp::coordinates(dat) = ~x +y #make sure its spatial

    #function to integrate ea_spline at each pixel
    mods = ithir::ea_spline(dat, var.name = label, vlow = 0, ...)
    f = function(x) integrate(approxfun(x$var.1cm, depths),
                            lower = min(dat), upper = max(dat), stop.on.error = F)$value

    #apply function to stack of spatrasters
    spR = terra::app(stack, f)
  }
  return(spR)
}
##END
