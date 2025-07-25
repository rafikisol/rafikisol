#'@title Add common soil and vegetative indices to Sentinel 2
#'
#'@name addIndices
#'
#'@author Trevan Flynn
#'
#'@description
#'Adds common soil and vegetative indices to sentinel 2A image collections. This should be conducted
#'after applying the s2Mask function to clean the images.
#'
#'@param s2Img Sentinel 2 ImageCollection to use.
#'
#'@return The Image Collection with the indices for each image in the series
#'
#'@export

addIndices = function(s2Img){

  #get bands and indices
  f1 = function(image){
    blue = image$select("B2")$rename("Blue")
    green = image$select("B3")$rename("Green")
    red = image$select("B4")$rename("Red")
    nir = image$select("B8")$rename("NIR")
    swir = image$select("B11")$rename("SWIR")
    ndvi = image$normalizedDifference(c('B8', 'B4'))$rename('NDVI')
    ndwi = image$normalizedDifference(c('B11', 'B8'))$rename('NDVI')
    bi = (red$pow(2)$multiply(blue$pow(2))$multiply(green$pow(2)))$
      divide(ee$Number(3)$sqrt())$rename("BI")
    ci = green$subtract(red)$divide(green$add(red))$rename('CI')
    ri = (red$pow(2))$divide(blue$multiply(green)$pow(3))$rename("RI")
    si = red$subtract(blue)$divide(red$add(blue))$rename('SI')

    #add bands to each image in the imageCollection
    image = image$addBands(blue)$addBands(green)$addBands(red)$addBands(nir)$
      addBands(swir)$addBands(ndvi)$addBands(ndwi)$addBands(bi)$addBands(ci)$
      addBands(ri)$addBands(si)
  }
  #apply function
  return(s2Img$map(f1))
}
