
#' rna_heatmap
#' 
#' @param mat matrix of values to be plotted as heatmap
#' @param x x xaxis labels, by default colnames of mat
#' @param y y axis labels, by default rownames of mat
#' @param row_order how to order rows?  See Details
#' @param col_order how to order columns? See Details
#' @param row_groups groups for rows, overriden if row_k is set and row_order is kmeans or hclust
#' @param col_groups groups for columns, overriden if col_k is set and col_order is kmeans or hclust
#' @param row_k number of clusters for rows, needed if row_order is kmeans or optional if hclust
#' @param col_k number of clusters for columns, needed if row_order is kmeans or optional if hclust
#' @param row_clust_dist distance function to use for row clustering if hierarchical clustering
#' @param col_clust_dist distance function to use for column clustering if hierarchical clustering
#' @param name Name for colorbar
#' @param source source name, useful for shiny
#' @param scale scale matrix by rows, cols or none
#' @param scale_method what method to use for scaling, either standardize, center, normalize
#' @param x_labels axis labels for x axis (default is x)
#' @param y_labels axis labels for y axis (default is NULL)
#' @param x_title x axis title
#' @param y_title y axis title (default is "Genes")
#' @param ... additional argument to make_main_hm
#' @return iHeatmap object, which can be printed or passed to \code{\link{plot_iHeatmap}} to
#' generate an interactive graphic
#' @export
#' @author Alicia Schep
rna_heatmap <- function(mat, 
                           x = iHeatmapR:::default_x(mat),
                           y = iHeatmapR:::default_y(mat),                   
                           row_order = c("hclust","none","kmeans","groups"),
                           col_order = c("groups","none","hclust","kmeans"),
                           row_groups = NULL,
                           col_groups = NULL,
                           row_k = NULL,
                           col_k = NULL,
                           row_clust_dist = stats::dist,
                           col_clust_dist = stats::dist,
                           name = "RNA-seq",
                           source = "HM",
                           scale = c("rows","cols","none"),
                           scale_method = c("standardize","center","normalize"),
                           x_labels = NULL,
                           y_labels = y,
                           x_title = NULL,
                           y_title = "Genes",
                           ...){
  
  row_order = match.arg(row_order)
  col_order = match.arg(col_order)
  scale = match.arg(scale)
  scale_method = match.arg(scale_method)
  
  iHeatmapR::simple_heatmap(mat, 
                            x = x,
                            y = y,                   
                            row_order = row_order,
                            col_order = col_order,
                            row_groups = row_groups,
                            col_groups = col_groups,
                            row_k = row_k,
                            col_k = col_k,
                            row_clust_dist = row_clust_dist,
                            col_clust_dist = col_clust_dist,
                            name = name,
                            source = source,
                            scale = scale,
                            scale_method = scale_method,
                            x_labels = x_labels,
                            y_labels = y_labels,
                            x_title = x_title,
                            y_title = y_title,
                            ...)
  
}

#' add_rna_heatmap
#' 
#' @param mat matrix of values to be plotted as heatmap
#' @param x x xaxis labels, by default colnames of mat
#' @param col_order how to order columns? See Details
#' @param col_groups groups for columns, overriden if col_k is set and col_order is kmeans or hclust
#' @param col_k number of clusters for columns, needed if row_order is kmeans or optional if hclust
#' @param col_clust_dist distance function to use for column clustering if hierarchical clustering
#' @param name Name for colorbar
#' @param scale scale matrix by rows, cols or none
#' @param scale_method what method to use for scaling, either standardize or center
#' @param x_labels axis labels for x axis (default is x)
#' @param y_labels axis labels for y axis (default is NULL)
#' @param x_title x axis title
#' @param y_title y axis title
#' @param buffer amount of space to leave empty before this plot, relative to size 
#' of first heatmap
#' @export
#' @author Alicia Schep
add_rna_heatmap <- function(p,
                               mat, 
                               x = iHeatmapR:::default_x(mat),
                               col_order = c("groups","none","hclust","kmeans"),
                               col_groups = NULL,
                               col_k = NULL,
                               col_clust_dist = stats::dist,
                               name = "Signal",
                               scale = c("rows","cols","none"),
                               scale_method = c("standardize","center","normalize"),
                               x_labels = x,
                               y_labels = NULL,
                               x_title = NULL,
                               y_title = NULL,
                               buffer = 0.2,
                               ...){
  
  
  col_order = match.arg(col_order)
  scale = match.arg(scale)
  scale_method = match.arg(scale_method)
  
  p %>% iHeatmapR::add_simple_heatmap(mat, 
                                  x = x,
                                  row_order = row_order,
                                  col_order = col_order,
                                  col_groups = col_groups,
                                  col_k = col_k,
                                  col_clust_dist = col_clust_dist,
                                  name = name,
                                  scale = scale,
                                  scale_method = scale_method,
                                  x_labels = x_labels,
                                  x_title = x_title,
                                  ...)
  
}