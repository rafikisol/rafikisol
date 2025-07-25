#'@title Function to calculate probabilities
#'
#'@name probCalc
#'
#'@author Trevan Flynn
#'
#'@description
#'Obtain the probabilities from rafiki() or rafikboot() algorithms. This function
#'also gets the data ready for a confusion index.
#'
#'@param models The results from rafiki() or rafikBoot().
#'@param covar The covariates used in the boostrapped models.
#'@param levels The number of classes being predicted.
#'
#'@return A list of probabilities for each bootstrap..
#'@export


#probabiliy predictions function
probCalc = function(models, covar, levels){

  #set models to probability
  pMods = lapply(models[["models"]],
                 function(x) {x$setOutputMode("MULTIPROBABILITY")})

  #create list of probs
  pPreds = lapply(pMods,
                  function(x) {covar$classify(x, "probability")$slice(0)})

  #make list to store image of each class (lists of lists)
  probs = list()

  #loop around classes
  for(i in 1:levels){
    #apply function
    probs[[i]] = lapply(pPreds, function(x) x$toArray(1)$arrayGet(c(i,0)))
  }

  return(probs)
}

