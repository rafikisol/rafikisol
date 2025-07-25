#'@title Calculate a confusion index (CI)
#'
#'@name ciCalc
#'
#'@author Trevan Flynn
#'
#'@description
#'Calculates the spatial CI from a bootstrapped ee model from Burrough et al. (1997)
#'where CI = 1 - (Pm - pm-1).Pm is the highest probability where Pm-1 is the second highest probability.
#'The final CI is the median of both Images. This function calculates the median values.
#'
#'@param probs Is the output from the probCalc() function.
#'
#'@returns An image of the confusion index.
#'
#'@references Burrough, P.A., Van Gaans, P.F.M., Hootsmans, R., (1997). Continuous classification in soil survey: spatial correlation, confusion and boundaries. Geoderma 77, 115–135.
#'https://doi.org/10.1016/S0016-7061(97)00018-9.
#'
#'@export
ciCalc = function(probs){

  #sort
  sortClss = lapply(probs, function(x) ee$ImageCollection(x)$
                      sort("array", FALSE))

  #first probs
  firstProb = ee$ImageCollection(lapply(sortClss, function(x)
    ee$ImageCollection(x)$max()))$ #take max of all
    median() #make image by taking the median

  #apply function over lists and take median
  secondProb = ee$ImageCollection(lapply(sortClss, function(x){
    x$map(function(image){
      return(image$updateMask(image$lt(firstProb)))})$
      max()}))$median()


  #calculate CI
  ci = ee$Image$constant(1)$subtract(firstProb$subtract(secondProb))

  #return CI
  return(ci$rename("CI"))
}

