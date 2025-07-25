#'@title Bootstapped Google Earth Engine classifiers
#'
#'@name rafikBoot
#'
#'@author Trevan Flynn
#'
#'@description
#'Bootstrapped Google Earth Engine classification models with either predifined
#'aggregation functions, custom functions or raw predictions.
#'
#'@param x Is the FeatureCollection with the sampling observations and labels.
#'@param y Is the Image of the covariates, where each band is a covariate.
#'@param label A character of what is to be predicted.
#'@param model Character of which model to use.
#'@param mode Should it be a classification, multiprobability, binomial probability or regression
#'@param fun Which aggregation function to use.
#'@param boot Integer of how many bootstrap iterations to perform.
#'@param fact The fraction of subsamples to perform.
#'@param scale What scale to perform the calculations in.
#'@param ... Additional arguments for the algorithm parameters (e.g., numberOfTrees)
#'
#'@export

rafikBoot = function(x, y, label, model, mode = "classification", fun = "mode", boot, fact, scale = y$projection()$nominalScale(), ...){

  #define which model
  if(model == "nb"){
    model = ee$Classifier$smileNaiveBayes
  }
  if(model == "gtb"){
    model = ee$Classifier$smileGradientTreeBoost
  }
  if(model == "svm"){
    model = ee$Classifier$libsvm
  }
  if(model == "cart"){
    model = ee$Classifier$smileCart
  }
  if(model == "maxent"){
    model = ee$Classifier$amnhMaxent
  }
  else{model = ee$Classifier$smileRandomForest
  }

  #make lists
  samps = list()

  #bootstrap samples
  for(i in 1:boot){
    samps[[i]] = x$randomColumn(seed = i)
  }

  #get training data
  training = lapply(samps, function(samps) samps$filter(ee$Filter$lte('random', fact)))

  #get covariates
  reSamps = lapply(training, function(training) y$sampleRegions(
    collection = training,
    properties = list(label),
    scale = scale,
    geometries = T))

  #run model
  if(mode == "classification"){
    mod = lapply(reSamps, function(reSamps) {model(...)$
        train(features = reSamps,
              classProperty = label,
              inputProperties = y$bandNames())})
  }

  if(mode == "probabilities"){
    mod = lapply(reSamps, function(reSamps) {model(...)$
        setOutputMode("MULTIPROBABILITY")$
        train(features = reSamps,
              classProperty = label,
              inputProperties = y$bandNames())})
  }

  if(mode == "binomial"){
    mod = lapply(reSamps, function(reSamps) {model(...)$
        setOutputMode("PROBABILITY")$
        train(features = reSamps,
              classProperty = label,
              inputProperties = y$bandNames())})
  }

  if(mode == "regression"){
    mod = lapply(reSamps, function(reSamps) {model(...)$
        setOutputMode("REGRESSION")$
        train(features = reSamps,
              classProperty = label,
              inputProperties = y$bandNames())})
  }

  #predict on each bootstrap
  pred = lapply(mod, function(mod) y$clip(region)$classify(mod, "classification"))

  #if maxent pull the probability density
  if(model == "maxent"){
    pred = lapply(pred, function(pred) pred$select(0))
  }

  #return object of function
  if(fun == "raw"){
    results = list(pred$rename("raw"), mod, reSamps)%>%
      setNames(c("predictions", "models", "samples"))

    return(results)
  }
  if(fun == "mode"){
    pred = ee$ImageCollection(pred)$
      mode()$
      rename("mode")
    return(list(pred, mod, reSamps)%>%
             setNames(c("predictions", "models", "samples")))
  }
  if(fun == "count"){
    pred = ee$ImageCollection(pred)$
      count()$
      rename("count")
    return(list(pred, mod, reSamps)%>%
             setNames(c("predictions", "models", "samples")))
  }
  if(fun == "mean"){
    pred = ee$ImageCollection(pred)$
      mean()$
      rename("mean")
    return(list(pred, mod, reSamps)%>%
             setNames(c("predictions", "models", "samples")))
  }
  if(fun == "median"){
    pred = ee$ImageCollection(pred)$
      median()$
      rename("median")
    return(list(pred, mod, reSamps)%>%
             setNames(c("predictions", "models", "samples")))
  }
  if(fun == "sum"){
    pred = ee$ImageCollection(pred)$
      sum()$
      rename("sum")
    return(list(pred, mod, reSamps)%>%
             setNames(c("predictions", "models", "samples")))
  }
  else(
    pred = fun(ee$ImageCollection(pred))$
      rename("predictions"))
  return(list(pred, mod, reSamps)%>%
           setNames(c("predictions", "models", "samples")))
}
