#'@title Oblique geographic coordinates (OGC) through GEE
#'
#'@name ogcCalc
#'
#'@author Trevan Flynn
#'
#'@description
#'An algorithm to calculate OGC in GEE originally developed by Møller et al. (2020).
#'This algorithm is largly taken from them, however, this algorithm does not use any
#'raster image as an input.It also allows to adjust the resolution and coordinate system accordingly.
#'We have found coarse coordinates as a covarate gives artifacts and smoothing helps this.
#'@param area A Feature, FeatureCollection or Geometry of the area of interest.
#'@param crs The coordinate system to use (best in projected) but not required.
#'@param scale The resolution to use.
#'@param nr The number of angles to produce.
#'
#'@return An image with each band representing the angles calculated
#'
#'@references Møller, A. B., Beucher, A. M., and Pouladi, N., Greve, M. H. (2020) Oblique geographic  coordinates as covariates for digital soil mapping, SOIL 6(2) 269 - 289https:/10.5194/soil-6-269-2020
#'
#'@export

ogcCalc = function(area, crs = "EPSG: 3857", scale = 250, nr = 6){

  #create x and y raster
  img = ee$Image$pixelCoordinates(crs)$ #x and y from coordinates
    clip(area)$ #clip to region
    reproject(crs =crs, scale = scale)

  #select x and y
  x = img$select("x")
  y = img$select("y")

  #angles and position
  pis <- pi*seq(0, 1, 1/nr)[2:(nr)] #sequence of angles
  pos <- seq(2, nr, 1)

  #need to edit locations for even numbers
  is.even <-function(x, tol = .Machine$double.eps^0.5){
    abs(x - round(x)) < tol
  }

  if(is.even(nr/2)) {
    pis <- pis[ -(nr/2)]
    pos <- pos[ -(nr/2)]
  }

  #Make a list for each image calculation
  calc1 = list()
  calc2 = list()
  ogc = list() #list of ogc

  #calculation 1
  for(i in 1:length(pis)){

    #get constant image of angles
    calc1[[pos[i]]] = ee$Image$constant(ee$Number(pis[i]))$
      clip(region)

    #calc 1
    calc1[[pos[i]]] = calc1[[pos[i]]]$subtract((y$divide(x))$atan())$cos()

    #calc 2
    calc2 = (x$pow(2)$add(y$pow(2)))$sqrt()

    #calc 3
    ogc[[pos[i]]] = calc1[[pos[i]]]$multiply(calc2)
  }

  #We put x and y back in because they can be useful.
  ogc[[1]] = x
  ogc[[(nr/2)+1]] = y

  #create labels (cant have "." in GEE names so substitute out for "_")
  labs=  paste0("pi", format(seq(0, 1, 1/nr)[1:nr], digits = 1, nsmall = 2))%>%
    gsub("[.]", "_",.)

  #add xy back in (otherwise names will be wrong)
  ogc = ee$ImageCollection(ogc)$
    toBands()$
    rename(labs)

  #return OGC
  return(ogc)
}
