#'@title Cloud and shadow mask for Sentinel 2A surface reflectance images from Google Earth Engine.
#'
#'@name s2Mask
#'
#'@author Trevan Flynn
#'
#'@description
#'Masks clouds and shadows as well as some angle corrections for each image of a Sentinel 2A image collection.
#'The function has no parameters.
#'
#'@param s2Img The Sentinel 2A SR ImageCollection to clean.
#'
#'@return An ImageCollection with clouds and shadows masked.
#'@export


s2Mask = function(s2Img){

  #Make function to map
  f1 = function(image){
    cloudProb = image$select('MSK_CLDPRB')
    snowProb = image$select('MSK_SNWPRB')
    cloud = cloudProb$lt(5)
    snow = snowProb$lt(5)
    scl = image$select('SCL')
    shadow = scl$eq(3)
    cirrus = scl$eq(10)
    mask = (cloud$And(snow))$And(cirrus$neq(1))$And(shadow$neq(1))
    return(image$updateMask(mask))
  }

  #map function
  return(s2Img$map(f1))
}
