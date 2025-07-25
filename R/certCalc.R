#'@title Calculate uncertainties from bootstrap ee models
#'
#'@author Trevan Flynn
#'
#'@description
#'Calculates the uncertainties from a bootstrap regression from rafikReg() or rafikBootReg().
#'
#'@param models Are the results from rafikReg() or rafikBootReg().
#'@param uncert Is the type of uncertainty one wants calculated. Either "range" (prediction range),
#'"stdev" (standard deviation), "var" (variance) or "se" (standard error).
#'@param limits A vector of confidence intervals (defaults to 5 and 95)
#'
#'@return An Image of the results
#'@export

certCalc = function(models, uncert = "range", limits = c(5,95)){

  if(uncert == "range"){

    #get prediction intervals
    interVals = ee$ImageCollection(models[[2]])$
      reduce(reducer = ee$Reducer$percentile(limits),
             parallelScale = 16)

    #calculate range
    predRange = interVals$select(1)$subtract(interVals$select(0))$
      rename("range")

    #Combine all three
    results = ee$ImageCollection(c(predRange, interVals))$
      toBands()$
      rename(c("range", "lowerLim", "upperLim"))

    return(results)
  }

  #variance reducer
  if(uncert == "var"){

    predVar = ee$ImageCollection(models[[2]])$
      reduce(reducer = ee$Reducer$variance(),
             parallelScale = 16)$
      rename("var")

    return(predVar)
  }

  #standard deviation reducer
  if(uncert == "stdev"){

    stdev = ee$ImageCollection(models[[2]])$
      reduce(reducer = ee$Reducer$stdDev(),
             parallelScale = 16)$
      rename("stdev")

    return(stdev)
  }

  if(uncert == "se")
    #standard dev
    stdev = ee$ImageCollection(models[[2]])$
      reduce(reducer = ee$Reducer$stdDev(),
             parallelScale = 16)$
      rename("stdev")

  #standard error
  se = stdev$multiply(ee$Number(qnorm(0.95)))$
    rename("se")

  return(se)
}
