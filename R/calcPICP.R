#'@title Prediction interval coverage probability
#'
#'@name calcPICP
#'
#'@author Trevan Flynn
#'
#'@description
#'A function to calculate the PICP of bootstrap ee models. It is meant to evaluate
#'uncertainties of the predictions. It is basically a complete replica of the code
#'by Malone et al., (2011), we just do most of the calculations on images first.
#'
#'@param models The results from regScale() or regBoot().
#'@param eval Validation data in as a Feature.
#'
#'@return Returns a list with a plot of PICP to confidence intervals, the results and
#'Lins concordance coefficient (how well the plot is in a 1:1 line).
#'
#'@references B.P. Malone, A.B. McBratney, B. Minasny, (2011). Empirical estimates of uncertainty for mapping continuous depth functions of soil attributes. Geoderma, Volume 160, Issues 3–4,Pages 614-626. https://doi.org/10.1016/j.geoderma.2010.11.013.
#'
#'@export

calcPICP = function(models, eval, label){

  #upper intervals
  qu = c(99.5,98.75,97.5, 95, 90, 80, 70, 60, 55, 52.5)

  #Calculate upper lims
  upperLims = ee$ImageCollection(models[[2]])$
    reduce(reducer = ee$Reducer$percentile(qu),
           parallelScale = 16)

  #lower intervals
  ql = c(2.5, 5, 10, 15, 20, 30, 45, 45.5, 46.75, 47.5)

  #calculate lower lims
  lowerLims = ee$ImageCollection(models[[2]])$
    reduce(reducer = ee$Reducer$percentile(ql),
           parallelScale = 16)

  #sample upper
  upS = upperLims$
    sampleRegions(
      collection = eval,
      properties = list(label),
      scale = models[[1]]$projection()$nominalScale(),
      tileScale = 16,
      geometries = F)%>%
    ee_as_sf(maxFeatures = 1e16)%>%
    sf::st_drop_geometry()

  #sample lower
  loS = lowerLims$
    sampleRegions(
      collection = eval,
      properties = list(label),
      scale = models[[1]]$projection()$nominalScale(),
      tileScale = 16,
      geometries = F)%>%
    ee_as_sf(maxFeatures = 1e16)%>%
    sf::st_drop_geometry()

  #convert to matrix and sort correctly
  samps = as.matrix(loS[, 1])#keep response
  loS = as.matrix(loS[, -c(1)])#remove response
  upS = as.matrix(upS[, -c(1)])

  #make matrix for intervals
  Mat<-matrix(NA, nrow=nrow(samps), ncol=length(qu))

  #loop through
  for(i in 1:ncol(Mat)){
    Mat[, i] = as.numeric(samps <= upS[, i] &
                            samps >= loS[, i])
  }

  #confidence levels
  cl<- c(99,97.5,95,90,80,60,40,20,10,5)

  #calculate PICP
  picp = (colSums(Mat)/nrow(Mat))*10

  #get data frame to plot
  results = data.frame(picp, cl)

  #see how it does with a straight line
  ccc = yardstick::ccc_vec(cl, picp)

  #make plot
  p = ggplot2::ggplot(data = results, ggplot2::aes(x = cl, y = picp))+
    ggplot2::geom_point(color = 'black')+
    ggplot2::geom_smooth(se = F, color = "red", method = 'lm', linetype = "dotted")+
    ggplot2::geom_abline(intercept = 0, slope = 1, color = 'blue')+
    ggplot2::labs(x = "Confidence level", y = "PICP")+
    ggplot2::xlim(c(0,1)) + ggplot2::ylim(c(0, 1))+
    ggplot2::theme_bw()

 return(setNames(list(p, ccc, results),
                  c("plot", "CCC", "results")))
}
