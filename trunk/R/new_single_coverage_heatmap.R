#' new_single_coverage_heatmap
#' @param mat coverage matrix
#' @param x x axis labels
#' @param y y axis labels
#' @param row_order row order method
#' @param k k to use for kmeans clustering or cutting heirarchical clustering
#' @param groups pre-determined groups for rows
#' @param signal signal along row
#' @param name name of colorbar
#' @param summary make summary plot, boolean
#' @param source source name in plotly
#' @return list with three components, 1) plot -- function to plot
#' 2) row_order -- row order used, 3) dendro hclust object if heirarchical clustering
#' performed
#' @export
#' @author Alicia Schep
new_single_coverage_heatmap <- function(mat, 
                                    x = iHeatmapR:::default_x(mat),
                                    y = iHeatmapR:::default_y(mat),                    
                                    row_order = c("none","hclust","kmeans","groups","signal"),
                                    k = NULL,
                                    groups = NULL,
                                    clust_dist = stats::dist,
                                    signal = NULL,
                                    name = "Signal",
                                    summary = TRUE,
                                    source = "HM",
                                    scale = c("rows","none"),
                                    scale_method = c("normalize","standardize","center"),
                                    ...){
  
  # TO DO: Add argument check
  
  row_order = match.arg(row_order)
  scale <- match.arg(scale)
  scale_method <- match.arg(scale_method)
  
  if (scale != "none") mat <- iHeatmapR:::scale_mat(mat, scale = scale, scale_method = scale_method)
  
  dendro = NULL
  
  if (row_order == "hclust"){
    dendro = flashClust::hclust(clust_dist(mat))
    row_order = dendro$order
    if (!is.null(k)){
      groups = cutree(dendro, k = k)
    }
  } else if (row_order == "kmeans"){
    stopifnot(!is.null(k))
    groups = kmeans(mat, centers = k)$cluster
    row_order = order(groups)
  } else if (row_order == "groups"){
    row_order = order(groups)
  } else if (row_order == "signal"){
    stopifnot(!is.null(signal))
    row_order = order(signal, decreasing = FALSE)
  } else{
    row_order = 1:nrow(mat)
  }
  
  p <- make_main_hm(mat, 
                    x = x, 
                    y = y, 
                    row_order = row_order, 
                    name = name,
                    y_pos = "none",
                    x_pos = "none",
                    ...) 
  
  if (!is.null(groups)){
    p <- p %>% add_row_groups(groups, side = "left")
  }
  if (!is.null(dendro)){
    p <- p %>% add_row_dendro(dendro, side = "left")
  }    
  if (!is.null(signal)){
    p <- p %>% add_row_signal(signal, "Total Signal")
  }
  if (summary){
    p <- p %>% add_col_summary(groups, showlegend = FALSE)
  }
  
  p$row_groups <- groups
  
  return(p)
}

#' add_coverage_heatmap
#' @param mat coverage matrix
#' @param x x axis labels
#' @param groups pre-determined groups for rows
#' @param signal signal along row
#' @param name name of colorbar
#' @param summary make summary plot, boolean
#' @param scale "rows" or "none"
#' @param scale_method "normalize", "standardize", "center"
#' @return iHeatmap object
#' @export
#' @author Alicia Schep
add_coverage_heatmap <- function(p,
                                 mat, 
                                 x = iHeatmapR:::default_x(mat),                    
                                 groups = p$row_groups,
                                 signal = NULL,
                                 name = "Signal",
                                 summary = TRUE,
                                 scale = c("rows","none"),
                                 scale_method = c("normalize","standardize","center"),
                                 ...){
  
  
  scale <- match.arg(scale)
  scale_method <- match.arg(scale_method)
  
  if (scale != "none") mat <- iHeatmapR:::scale_mat(mat, scale = scale, scale_method = scale_method)
  
  dendro = NULL
  
  p <- p %>% add_main_hm(mat, 
                    x = x, 
                    name = name,
                    y_pos = "none",
                    x_pos = "none",
                    ...) 
  
  if (!is.null(signal)){
    p <- p %>% add_row_signal(signal, "Total Signal")
  }
  if (summary){
    p <- p %>% add_col_summary(groups, showlegend = FALSE)
  }
  
  return(p)

}


#' multi_coverage_heatmap
#' @param mat coverage matrix
#' @param x x axis labels
#' @param y y axis labels
#' @param row_order row order method
#' @param k k to use for kmeans clustering or cutting heirarchical clustering
#' @param groups pre-determined groups for rows
#' @param signal signal along row
#' @param name name of colorbar
#' @param summary make summary plot, boolean
#' @param source source name in plotly
#' @param scale scale rows? or none
#' @param scale_method method to use for scaling
#' @param share_z share colorbar between heatmaps
#' @return iHeatmap object
#' @export
#' @import iHeatmapR
#' @author Alicia Schep
multi_coverage_heatmap <- function(mats, 
                                        x = iHeatmapR:::default_x(mats[[1]]),
                                        y = iHeatmapR:::default_y(mats[[1]]),                    
                                        row_order = c("none","hclust","kmeans","groups","signal"),
                                        k = NULL,
                                        groups = NULL,
                                        clust_dist = stats::dist,
                                        clust_include = 1, 
                                        signal = NULL,
                                        name = "Signal",
                                        summary = TRUE,
                                        source = "HM",
                                        scale = c("rows","none"),
                                        scale_method = c("normalize","standardize","center"),
                                   share_z = TRUE,
                                   
                                        ...){
  
  # TO DO: Add argument check
  
  row_order = match.arg(row_order)
  scale <- match.arg(scale)
  scale_method <- match.arg(scale_method)
  
  if (scale != "none"){
    mats <- lapply(mats, iHeatmapR:::scale_mat, scale, scale_method)
  } 
  
  if (share_z){
    zmax <- max(sapply(mats, max, na.rm = TRUE))
    zmin <- min(sapply(mats, min, na.rm = TRUE))
  }
  
  dendro = NULL
  
  if (row_order == "hclust"){
    dendro = flashClust::hclust(clust_dist(do.call(cbind,mats[clust_include])))
    row_order = dendro$order
    if (!is.null(k)){
      groups = cutree(dendro, k = k)
    }
  } else if (row_order == "kmeans"){
    stopifnot(!is.null(k))
    groups = kmeans(do.call(cbind,mats[clust_include]), centers = k)$cluster
    row_order = order(groups)
  } else if (row_order == "groups"){
    row_order = order(groups)
  } else if (row_order == "signal"){
    stopifnot(!is.null(signal))
    row_order = order(signal, decreasing = FALSE)
  } else{
    row_order = 1:nrow(mats[[1]])
  }
  
  if (share_z){
    p <- make_main_hm(mats[[1]], 
                      x = x, 
                      y = y, 
                      row_order = row_order, 
                      name = name,
                      y_pos = "none",
                      x_pos = "none",
                      zmin = zmin,
                      zmax = zmax,
                      ...) 
  } else {
    p <- make_main_hm(mats[[1]], 
                    x = x, 
                    y = y, 
                    row_order = row_order, 
                    name = names(mats)[1],
                    y_pos = "none",
                    x_pos = "none",
                    ...)     
  }
  

  if (summary){
    p <- p %>% add_col_summary(groups, showlegend = FALSE)
  }  
  if (!is.null(signal)){
    p <- p %>% add_row_signal(signal, "Total Signal")
  }
  
  if (length(mats) > 1){
    for (i in 2:length(mats)){
      if (share_z){
        p <- p %>% add_main_hm(mats[[i]], 
                               x = x, 
                               y_pos = "none",
                               x_pos = "none",
                               show_colorbar = FALSE,
                               zmax = zmax,
                               zmin = zmin,
                               ...) 
      } else{
        p <- p %>% add_main_hm(mats[[i]], 
                             x = x, 
                             name = names(mats)[i],
                             y_pos = "none",
                             x_pos = "none",
                             ...) 
      }
      
      
      if (!is.null(signal)){
        p <- p %>% add_row_signal(signal, "Total Signal")
      }
      if (summary){
        p <- p %>% add_col_summary(groups, showlegend = FALSE, yaxis = "y2")
      }
      
    }
  }
  
  
  if (!is.null(groups)){
    p <- p %>% add_row_groups(groups, side = "left")
  }
  if (!is.null(dendro)){
    p <- p %>% add_row_dendro(dendro, side = "left")
  }    

  return(p)
}



add_dist_to_tss <- function(ranges,
                            tss,
                            ...){
  
  
}




