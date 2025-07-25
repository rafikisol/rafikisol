#'@title Evaluate regression models in GEE
#'
#'@name evalCalc

evalCalc = function(models, eval, label){

  pred = models[[1]]$
    sampleRegions(
      collection = eval,
      properties = list(label),
      tileScale = 16,
      geometries = F)$
    select(c(label, "median"))

  #Pearson's correlation coefficient
  r = pred$reduceColumns(
    reducer = ee$Reducer$pearsonsCorrelation(),
    selectors = list(label, "median"))

  #R2
  r2 = ee$Number(r$get("correlation"))$ #get correlation
  pow(2)$getInfo() #square it

  addResid <- function(feature) {
    res <- ee$Number(feature$get(label))$ #subtract observed from predicted
      subtract(ee$Number(feature$get("median")))
    feature$set(list(res = res)) #create new feature in featureCollection
  }

  #apply function to FeatureCollection for residuals
  res = pred$map(addResid)

  #calculate RMSE
  rmse = ee$Array(res$aggregate_array("res"))$
    pow(2)$ #square it
    reduce("mean", list(0))$ #get mean for the residuals
    sqrt()$getInfo()

  return(list(r2, rmse))
}
