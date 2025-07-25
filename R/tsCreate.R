#'@title Create a time series that can be plotted in R from an image collection
#'
#'@name tsCreate
#'
#'@author Trevan Flynn
#'
#'@description
#'Reduces the image collection into a feature which can be transformed into an sf object
#'to plot in R and for further analysis. Please note that geometry must be labeled geometry and that
#'it reduces the images based on the median.
#'
#'@param imgCol ImageCollection to reduce (e.g., Sentinel time series).
#'@param reducer The reducer to be used. Either a character or a built reducer (e.g., ee$Reducer$mean()).
#'@param geom Point or polygon of location to be summarised.
#'@param scale Scale at which to reduce values. Defaults to image$projection()$nominalScale().
#'@param dateFormat Date format you want the dataframe to come out as.
#'@param sf Logical, should the data frame be converted to an sf object into local space.
#'
#'@return A feature with values of image for each date in time series.
#'@export

tsCreate = function(imgCol, reducer = "mean", scale, geom,
                    dateFormat = "YYYY-MM-dd", sf = FALSE){

  #Make a function to map over image collection
  f1 = function(image){

    value = image$reduceRegion(reducer = reducer, geometry = geom, scale = scale)
    date = ee$List(image$date()$format(dateFormat))
    ft = ee$Feature(NULL, value)
    ft = ft$set("Time", date)

    return(ft)
  }

  #apply function on image collection
  fc = ee$FeatureCollection(imgCol$map(f1))

  if(sf == TRUE){
    return(ee_as_sf(fc))
  }

  return(fc)
}
#
