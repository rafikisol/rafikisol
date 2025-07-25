#'@title Classify and evaluate catSpline
#'@name evalSpline
#'@author Trevan Flynn
#'
#'@description Takes the output of catSpline() and puts out categorical evaluation
#'statistics from the observed layers to the fitted spline.
#'
#'@param spl is the spline output.
#'
#'@returns a list from the raw probabilities but classified
#'@export

evalSpline = function(spl){
  lookup <- spl$lookup
  fit <- spl$obs.preds

  #Confusion matrix
  conMat <- table(Obvserved = fit$Obvs, Predicted = fit$Pred, useNA = "ifany")

  #Kappa
  Po <- sum(diag(conMat))/sum(conMat)
  rtot <- rowSums(conMat)
  ctot <- colSums(conMat)
  Pe<- sum((rtot*ctot)/sum(conMat)^2)

  kappa <- round((Po - Pe) / (1 - Pe), 2)

  #accuracies
  OA <- ceiling(sum(diag(conMat)) / sum(conMat)* 100)
  PA <- ceiling(diag(conMat)/colSums(conMat) * 100)
  UA <- ceiling(diag(conMat)/rowSums(conMat) * 100)

  #Gather results
  reSlt = list(OA, kappa, PA, UA, conMat, fit)
  names(reSlt) = c("Overall_accuracy", "Kappa", "Producers_accuracy", "Users_accuracy", "Confusion_matrix", "Raw_fit")
  return(reSlt)
}
##END
