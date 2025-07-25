#'@title Evaluation stats for categorical data
#'
#'@name catEval
#'
#'@author Trevan Flynn
#'
#'@description
#'Evaluate categorical results from a bootstrapped ee model through rgee functions
#'If you would like more options look into the yardstick package.
#'
#'@param models The model results from catScale() or catBoot().
#'@param eval A FeatureCollection to evaluate the model.
#'@param response The response of the model (e.g., "texture").
#'
#'@return Returns a list containing the kappa coefficient, accuracy, consumer's
#'accuracy (CA) and producer's accuracy (PA).
#'@export

catEval = function(models, eval, response){

  #sample data
  eval = preds[['modal']]$sampleRegions(
    collection = eval,
    properties = list(response),
    scale = models[[1]]$projection()$nominalScale(),
    geometries = T)

  #confusion matrix
  errMat = eval$errorMatrix(response, "modal")

  #some metrics
  kap = errMat$kappa()$getInfo()%>% tidyr::as_tibble()
  acc = errMat$accuracy()$getInfo()%>% tidyr::as_tibble()
  ca = errMat$consumersAccuracy()$getInfo()%>% tidyr::as_tibble()
  pa = errMat$producersAccuracy()$getInfo()%>% tidyr::as_tibble()

  #list
  met = list(kap, acc, ca, pa)

  return(setNames(mets, c("kappa", "accuracy", "CA", "PA")))
}
