#'@title Function to bootstrap GEE classifiers
#'
#'@name catBoot
#'
#'@author Trevan Flynn
#'
#'@description
#'A function to bootstrap ee models which is useful for small datasets and gives
#'a more robust measure of uncertainty. It also gives a better estimate of variable
#'importance.
#'
#'@param covar The covariate Image.
#'@param response The FeatureCollection with what needs to be predicted.
#'@param label Character of the response to predict.
#'@param model Character of which ee model to use.
#'@param boot Integer of how many bootstraps.
#'@param ratio The ration of subsampling to use in each bootstrap.
#'
#'@returns A list of the modal predictions, all preds, all models and the training and
#'evaluation samples for each bootstrap.
#'@export

catBoot = function(covar, response, label, model = NULL, boot = 10, ratio = 0.6, ...){

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
  else{model = ee$Classifier$smileRandomForest
  }

  #extract all samples
  samps = covar$sampleRegions(
    collection = response,
    properties = list(label),
    scale = scale,
    geometries = T
  )

  #bootstrap the samples
  resSamp = list()
  training = list()
  testing = list()

  for(i in 1:boot){

    #add a random column to split each time
    ResSamps[[i]] = samps$randomColumn(seed = i, distribution = "uniform")

    #get training and test sets
    training[[i]] = samps[[i]]$filter(ee$Filter$lte('random',ratio))
    validation[[i]] = samps[[i]]$filter(ee$Filter$gt('random',ratio))
  }

  #build and train models
  mods = lapply(samps, function(x)
    {model(...)$train(features = x, classProperty = label,
                                                       inputProperties = covar$bandNames())})
  #predict models over covariate image
  preds = lapply(mods, function(x) {covar$classify(x, "predictions")})

  #get mode of the imageCollection
  predModal = ee$ImageCollection(preds)$
    mode()$
    rename("modal")

  #names list for ease of use
  return(setNames(list(predModal, preds, mods, training, testing),
                  c("modal","predictions", "models", "training", "testing")))

  }
