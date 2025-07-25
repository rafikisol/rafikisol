#'@title Calculate importance of bootstrapped ee models
#'
#'@name impCalc
#'
#'@author Trevan Flynn
#'
#'@description
#' Code to get the average variable importance from a bootstrapped model in GEE.
#'
#'@param models Is the results from a bootstrapped ee model (regScale(), catScale(), regBoot() or catBoot()).
#'@param mn The number of covariates to display in the plot (all others will be in the returned value).
#'@return A list of [[1]] a ggplot2 lollipop plot and [[2]] results of importance in a tibble.
#'
#'@export

impCalc = function(models, mn){

  if(is.list(models)){
  #extract models
  imp = lapply(models[["models"]],
               function(x){x$explain()$get("importance")$getInfo()%>%#get info
                   dplyr::as_tibble()%>% #make tibble
                   tidyr::pivot_longer(cols = c(1:ncol(.)))})%>%
    dplyr::bind_rows()%>% #bind rows
    dplyr::group_by(name)%>%
    dplyr::mutate(meanImp = mean(value))%>% #calculate by group
    dplyr::distinct(name, meanImp)%>% #only distinct
    dplyr::arrange(desc(meanImp))

  #add a ggplot of results
  p = ggplot2::ggplot(imp[1:mn, ], #get top rows
                      ggplot2::aes(y = meanImp, x= reorder(name, meanImp)))+
    ggplot2::geom_segment(ggplot2::aes(xend = name, yend = meanImp, y = 0, color = name),show.legend = F)+
    ggplot2::geom_point(ggplot2::aes(color = name), show.legend = F)+
    ggplot2::scale_color_viridis_d(option = "turbo")+
    ggplot2::labs(x = "Covariate", y = "Importance")+
    ggplot2::theme_bw()+
    ggplot2::coord_flip()
  }

  else{
    models$explain()$get("importance")$getInfo()%>%
      dplyr::as_tibble()%>%
      tidyr::pivot_longer(cols = c(1:ncol(.)))%>%
      dplyr::bind_cols()%>%
      dplyr::group_by(name)%>%
      dplyr::mutate(meanImp = mean(value))%>%
      dplyr::distinct(name, meanImp)%>%
      dplyr::arrange(desc(meanImp))

    p = ggplot2::ggplot(imp[1:mn, ], #get top rows
                        ggplot2::aes(y = meanImp, x= reorder(name, meanImp)))+
      ggplot2::geom_segment(ggplot2::aes(xend = name, yend = meanImp, y = 0, color = name),show.legend = F)+
      ggplot2::geom_point(ggplot2::aes(color = name), show.legend = F)+
      ggplot2::scale_color_viridis_d(option = "turbo")+
      ggplot2::labs(x = "Covariate", y = "Importance")+
      ggplot2::theme_bw()+
      ggplot2::coord_flip()
  }
  #return and set names
  return(setNames(list(p, imp), c("plot", "importance")))
}

