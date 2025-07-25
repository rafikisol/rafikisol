#'@title Confusion index (CI) for ee Models
#'
#'@name rafikCI
#'
#'@author Trevan Flynn
#'
#'@description
#'Calculate a confusion index for an ee$Classifier through Burrough et al. 1997.
#'Unlike the ciCalc function this calculates a CI for a single model.
#'
#'@param probs Probabilities from running setOutputMode("MULTIPROBABILITY") to
#'an ee$Classifier.
#'
#'@return A CI image.
#'
rafikCI = function(probs){

  #first prob
  firstProb = probs$max()

  #second prob - need to mask first
  secondProb = probs$map(function(image){
    return(image$updateMask(image$lt(firstProb)))
  })$max()

  #calculate CI
  ci = ee$Image$constant(1)$subtract(firstProb$subtract(secondProb))

  return(ci$rename("CI"))
}
