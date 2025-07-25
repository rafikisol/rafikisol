#' @title Downscale categorical maps through GEE
#'
#' @name rafikisol
#'
#' @author Trevan Flynn
#'
#' @description
#' A function to downscale categorical soil maps with ee models through Monte Carlo sampling.
#' The algorithm is basically the same as the dissever algorithm by
#' Malone et al. (2012) and Roudier et al. (2017) besides the fact it is for categorical data
#' and utilizes GEE and therefore, enables parallel processing and lowers computational power.
#'
#' @param covar The covariates to be used in the model (should be fine resolution).
#' @param response The image with the data to be sampled on (should be coarse resolution).
#' @param model A character of which ee model to be used being random forest ("rf"),
#' gradient tree boost ("gtb"), cart ("cart") or support vector machines ("svm").
#' @param stype Character of sampling type to perform either stratified ("s") or random ("r").
#' @param scale The scale at which to perform analysis. Defaults to scale of covar.
#' @param sampNum The number of samples for each bootstrap. For "s", its the amount per strata.
#' @param boot An integer of the number of bootstraps to perform.
#' @param ... Parameters to add to the chosen ee model (e.g., numberOfTrees).
#'
#' @references Brendan P. Malone, Alex B. McBratney, Budiman Minasny, Ichsani Wheeler, (2012). A general method for downscaling earth resource information. Computers & Geosciences, Volume 41, Pages 119-125. https://doi.org/10.1016/j.cageo.2011.08.021.
#' @references P. Roudier, B.P. Malone, C.B. Hedley, B. Minasny, A.B. McBratney, (2017). Comparison of regression methods for spatial downscaling of soil organic carbon stocks maps. Computers and Electronics in Agriculture. Volume 142, Part A, Pages 91-100. https://doi.org/10.1016/j.compag.2017.08.021.
#'
#' @returns Returns a list of four with the mean predictions, all predictions,
#' all models and all samples of the bootstraps.
#'
#' @export

#function
rafikisol = function(covar, response, scale = covar$projection()$nominalScale(),
                    model = NULL, stype = 'r', sampNum = 100, boot = 10,...){

  #list (only for bootstrap)
  resSamp = list()

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

  #loop around samples
  for(i in 1:boot) {

    #stratified samples
    if(stype == "s"){
      #stratified samples
      resSamp[[i]] = response$ceil()$int()$stratifiedSample(
        numPoints = sampNum,
        classBand = label,
        scale = response$projection()$nominalScale(),
        region = region,
        geometries = T,
        seed = i)
    }

    #random sampling
    else{
      #sample large raster
      resSamp[[i]] = response$int()$sample(
        region = region,
        numPixels = sampNum,
        dropNulls = T,
        scale = response$projection()$nominalScale(),
        seed = i,
        geometries = T)
    }
  }

  #sample fine covariates
  samps = lapply(resSamp, function(x) {covar$sampleRegions(
    collection = x,
    properties = list(label),
    scale = scale,
    geometries = T)})

  #build and train models
  mods = lapply(samps, function(x) {model(...)$train(features = x,
                                                     classProperty = label,
                                                     inputProperties = covar$bandNames())})

  #predict models over covariate image
  preds = lapply(mods, function(x) {covar$classify(x, "predictions")})

  #get mode of the imageCollection
  predModal = ee$ImageCollection(preds)$
    mode()$
    rename("modal")

  #names list for ease of use
  return(setNames(list(predModal,preds, mods, samps),
                  c("modal","predictions", "models", "samples")))
}

