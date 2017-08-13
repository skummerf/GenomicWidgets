#' @export
log10rowMeans <- function(x, na.rm = TRUE, pseudo = 1){
  log10(rowMeans(x, na.rm = na.rm) + pseudo)
}

#' coverage_heatmap
#' 
#' @param data single coverage matrix or list of coverage matrices
#' @param x x axis labels
#' @param y y axis labels
#' @param row_order row order method
#' @param k k to use for kmeans clustering or cutting heirarchical clustering
#' @param groups pre-determined groups for rows
#' @param clust_dist distance function to use for clustering
#' @param signal signal along row
#' @param plot_signal add an annotation heatmap showing the average signal? 
#' default is TRUE
#' @param name name of colorbar
#' @param summary make summary plot, boolean, default is TRUE
#' @param scale_method how to scale matrix before displaying in heatmap
#' @param pct percentile to use if scale_method is "PercentileMax"
#' @param scale_factor scale_factor to use if scale_method is "scalar"
#' @param show_xlabels show xlabels?  default is TRUE
#' @param start label for start of x range
#' @param end label for end of x range
#' @param xlab x axis label
#' @param layout list of layout attributes
#' @return iheatmap object
#' @export
#' @author Alicia Schep
#' @rdname coverage_heatmap
#' @name coverage_heatmap
#' @aliases coverage_heatmap,matrix-method coverage_heatmap,list-method
#' add_coverage_heatmap,IheatmapHorizontal,matrix-method
#' @export
#' @import iheatmapr
setGeneric("coverage_heatmap", 
           function(data, ...) standardGeneric("coverage_heatmap"))

#' @rdname coverage_heatmap
#' @export
setGeneric("add_coverage_heatmap", 
           function(p, data, ...) standardGeneric("add_coverage_heatmap"))

# Method for bamfile or list of bamfiles (or bigwig files)
setMethod("coverage_heatmap", c(data = "character"),
          function(data, 
                   windows,
                   ...){
            
            # ADD HERE
          })

# Method for ScoreMatrix
setMethod("coverage_heatmap", c(data = "ScoreMatrix"),
          function(data, 
                   start, 
                   end, 
                   x = seq(start, end, length.out = ncol(data)),
                   y = default_y(data), 
                   ...){
            
            coverage_heatmap(as.matrix(data),
                             x = x, y = y,
                             start = start,
                             end = end, ...)
          })


# Method for ScoreMatrixList
setMethod("coverage_heatmap", c(data = "ScoreMatrixList"),
          function(data, 
                   start, 
                   end, 
                   x = seq(start, end, length.out = ncol(data)),
                   y = default_y(data[[1]]), 
                   ...){
            
            data_list <- lapply(data, as.matrix)
            coverage_heatmap(data_list,
                             x = x, y = y,
                             start = start,
                             end = end, ...)
          })



setMethod("coverage_heatmap", c(data = "matrix"),
          function(data, 
                   x = default_x(data),
                   y = default_y(data), 
                   row_order = c("signal","hclust","kmeans","groups","none"),
                   k = NULL,
                   groups = NULL,
                   clust_dist = stats::dist,
                   signal = log10rowMeans(data),
                   plot_signal = TRUE,
                   name = "Coverage",
                   signal_name = "Avg. (log10)",
                   summary = TRUE,
                   scale_method = c("localRms", 
                                    "localMean", 
                                    "localNonZeroMean", 
                                    "PercentileMax", 
                                    "scalar", 
                                    "none"),
                   pct = 0.95,
                   scale_factor = 1, 
                   show_xlabels = TRUE,
                   start = x[1],
                   end = default_end(x),     
                   col_title = "Position",
                   layout = list(font = list(size = 10)),
                   ...){
            
            row_order <- match.arg(row_order)
            scale_method <- match.arg(scale_method)
                  
            
            if (scale_method != "none"){
              normdata <- normalize_coverage_matrix(data, 
                                               method = scale_method, 
                                               pct = pct, 
                                               scalar = scale_factor)
            } else{
              normdata <- data
            }
              
            p <- main_heatmap(normdata,
                              name = name,
                              x = x, 
                              y = y,
                              x_categorical = FALSE,
                              layout = layout,
                              ...)
            
            if (!is.null(groups) && row_order != "groups"){
              p <- add_row_groups(p, groups, side = "left")
            }     
            
            if (row_order == "signal"){
              p <- reorder_rows(p, order(signal))
            } else if (row_order != "none"){
              p <- add_row_clustering(p, 
                                      method = row_order,
                                      k = k, 
                                      groups = groups,
                                      clust_dist = clust_dist)
            }
            
            if (show_xlabels){
              if ("0" %in% x || 0 %in% x){
                ticktext = c(start, "0",end)
                tickvals = as.numeric(c(x[1], 0, x[length(x)]))
              } else{
                ticktext = c(start, end)
                tickvals = as.numeric(c(x[1], x[length(x)]))
              }
              
              p <- add_col_labels(p,
                                  ticktext = ticktext, 
                                  tickvals = tickvals)
              
            } 
            
            if (is.character(col_title) && nchar(col_title) > 0){
              p <- add_col_title(p,
                                 col_title, 
                                 side = "bottom")
            }
            
            if (plot_signal){
              p <- add_row_signal(p,
                                  signal, signal_name)
            }
            if (summary){
              p <- add_col_summary(p,
                                   groups = groups, showlegend = FALSE)
            }
            
            return(p)
              
})

setMethod("add_coverage_heatmap", c(p = "IheatmapHorizontal", data = "matrix"),
          function(p,
                   data, 
                   x = default_x(data),
                   signal = log10rowMeans(data),
                   plot_signal = TRUE,
                   name = "Coverage",
                   signal_name = "Avg. (log10)",
                   summary = TRUE,
                   scale_method = c("localRms", 
                                    "localMean", 
                                    "localNonZeroMean", 
                                    "PercentileMax", 
                                    "scalar", 
                                    "none"),
                   pct = 0.95,
                   scale_factor = 1, 
                   show_xlabels = TRUE,
                   start = x[1],
                   end = default_end(x),     
                   col_title = "Position",
                   ...){
            
            scale_method <- match.arg(scale_method)
            
            if (scale_method != "none"){
              normdata <- normalize_coverage_matrix(data, 
                                                    method = scale_method, 
                                                    pct = pct, 
                                                    scalar = scale_factor)
            } else{
              normdata <- data
            }
            
            p <- add_main_heatmap(p,
                                  normdata,
                                  name = name,
                                  x = x, 
                                  x_categorical = FALSE,
                                  ...)
           
            if (show_xlabels){
              if ("0" %in% x || 0 %in% x){
                ticktext = c(start, "0",end)
                tickvals = as.numeric(c(x[1], 0, x[length(x)]))
              } else{
                ticktext = c(start, end)
                tickvals = as.numeric(c(x[1], x[length(x)]))
              }
              
              p <- add_col_labels(p,
                                  ticktext = ticktext, 
                                  tickvals = tickvals)
              
            } 
            
            if (is.character(col_title) && nchar(col_title) > 0){
              p <- add_col_title(p,
                                 col_title, 
                                 side = "bottom")
            }
            
            if (plot_signal){
              p <- add_row_signal(p,
                                  signal, signal_name)
            }
            if (summary){
              p <- suppressWarnings(add_col_summary(p,
                                                    groups = TRUE, 
                                                    showlegend = FALSE))
            }
            
            return(p)
            
          })          
          
setMethod("coverage_heatmap", c(data = "list"),
          function(data, 
                   x = default_x(data[[1]]),
                   y = default_y(data[[1]]), 
                   row_order = c("signal","hclust","kmeans","groups","none"),
                   k = NULL,
                   groups = NULL,
                   clust_dist = stats::dist,
                   cluster_by = c("first","all"), 
                   signal = lapply(data, log10rowMeans),
                   plot_signal = TRUE,
                   name = "Coverage",
                   signal_name = "Avg. (log10)",
                   summary = TRUE,
                   scale_method = c("localRms", 
                                    "localMean", 
                                    "localNonZeroMean", 
                                    "PercentileMax", 
                                    "scalar", 
                                    "none"),
                   pct = 0.95,
                   scale_factor = 1, 
                   show_xlabels = TRUE,
                   start = x[1],
                   end = default_end(x),     
                   col_title = "Position",
                   layout = list(font = list(size = 10)),
                   ...){
            
            stopifnot(all(sapply(data, inherits, "matrix")))
            if (length(unique(lapply(data, nrow))) > 1) 
              stop("All input matrices must be of same length")
            
            row_order <- match.arg(row_order)
            scale_method <- match.arg(scale_method)
            cluster_by <- match.arg(cluster_by)
            
            if (scale_method != "none"){
              normdata <- normalize_coverage_matrix(data, 
                                                    method = scale_method, 
                                                    pct = pct, 
                                                    scalar = scale_factor)
            } else{
              normdata <- data
            }
            
            p <- main_heatmap(normdata[[1]],
                              name = name,
                              x = x, 
                              y = y,
                              x_categorical = FALSE,
                              layout = layout,
                              ...)
            
            if (!is.null(groups) && row_order != "groups"){
              p <- add_row_groups(p, groups, side = "left")
            }     
            
            if (cluster_by == "first"){
              if (row_order == "signal"){
                p <- reorder_rows(p, order(signal[[1]]))
              } else if (row_order != "none"){
                p <- add_row_clustering(p, 
                                        method = row_order,
                                        k = k, 
                                        groups = groups,
                                        clust_dist = clust_dist)
              }
            } else{
              if (row_order == "signal"){
                agg_signal <- Reduce("+",
                                     lapply(signal, function(z) z / sum(z)))
                p <- reorder_rows(p, order(agg_signal))
              } else if (row_order == "hclust"){
                dendro <- fastcluster::hclust(clust_dist(do.call(cbind,
                                                                 normdata)))
                if (!is.null(k)){
                  p <- add_row_groups(p, cutree(dendro, k = k),
                                      name = "Row<br>Clusters")
                }
                p <- add_row_dendro(p, dendro)
              } else if (row_order == "kmeans"){
                p <- add_row_clusters(p, 
                                      kmeans(do.call(cbind, normdata), 
                                             centers = k)$cluster )
              } else if (row_order == "groups"){
                p <- add_row_clustering(p, 
                                        method = "groups",
                                        groups = groups)
              }
            }
            
            
            if (show_xlabels){
              if ("0" %in% x || 0 %in% x){
                ticktext = c(start, "0",end)
                tickvals = as.numeric(c(x[1], 0, x[length(x)]))
              } else{
                ticktext = c(start, end)
                tickvals = as.numeric(c(x[1], x[length(x)]))
              }
              
              p <- add_col_labels(p,
                                  ticktext = ticktext, 
                                  tickvals = tickvals)
              
            } 
            
            if (is.character(col_title) && nchar(col_title) > 0){
              p <- add_col_title(p,
                                 col_title, 
                                 side = "bottom")
            }
            
            if (plot_signal){
              p <- add_row_signal(p,
                                  signal[[1]], 
                                  signal_name)
            }
            if (summary){
              summary_yaxis = "summary"
              p <- add_col_summary(p,
                                   groups = groups, 
                                   showlegend = FALSE,
                                   yname = summary_yaxis)
            }
            
            if (!is.null(names(data)[1]))
              p <- add_col_title(p, names(data)[1], side = "top")
            
            if (length(data) > 1){
              for (i in 2:length(data)){
                p <- add_main_heatmap(p,
                                      normdata[[i]], 
                                        x = x, 
                                      name = name,
                                        x_categorical = FALSE,
                                        ...) 
                
                if (plot_signal){
                  p <- add_row_signal(p,
                                            signal[[i]], 
                                            name = signal_name)
                }
                if (summary){
                  p <-add_col_summary(p,
                                      groups, 
                                      showlegend = FALSE, 
                                      yname = summary_yaxis)
                }
                
                if (show_xlabels){
                  p <- add_col_labels(p,
                                            ticktext = ticktext, 
                                      tickvals = tickvals)
                }
                if (!is.null(col_title)){
                  p <- add_col_title(p, col_title)
                }
                if (!is.null(names(data)[i]))
                  p <- add_col_title(p, names(data)[i], side = "top")
              }
            }
            
            return(p)
            
          })          


setMethod("add_coverage_heatmap", c(p = "IheatmapHorizontal", data = "list"),
          function(p,
                   data, 
                   x = default_x(data[[1]]),
                   signal = lapply(data, log10rowMeans),
                   plot_signal = TRUE,
                   name = "Coverage",
                   signal_name = "Avg. (log10)",
                   summary = TRUE,
                   scale_method = c("localRms", 
                                    "localMean", 
                                    "localNonZeroMean", 
                                    "PercentileMax", 
                                    "scalar", 
                                    "none"),
                   pct = 0.95,
                   scale_factor = 1, 
                   show_xlabels = TRUE,
                   start = x[1],
                   end = default_end(x),     
                   col_title = "Position",
                   ...){
            
            stopifnot(all(sapply(data, inherits, "matrix")))
            if (length(unique(lapply(data, nrow))) > 1) 
              stop("All input matrices must be of same length")
            
            scale_method <- match.arg(scale_method)

            if (scale_method != "none"){
              normdata <- normalize_coverage_matrix(data, 
                                                    method = scale_method, 
                                                    pct = pct, 
                                                    scalar = scale_factor)
            } else{
              normdata <- data
            }
            
            summary_yaxis <- paste0("summary", length(plots(p)))
            if (show_xlabels){
              if ("0" %in% x || 0 %in% x){
                ticktext = c(start, "0",end)
                tickvals = as.numeric(c(x[1], 0, x[length(x)]))
              } else{
                ticktext = c(start, end)
                tickvals = as.numeric(c(x[1], x[length(x)]))
              }
            }
            
            for (i in seq_along(data)){
              p <- add_main_heatmap(p,
                                    normdata[[i]], 
                                    x = x, 
                                    name = name,
                                    x_categorical = FALSE,
                                    ...) 
              
              if (plot_signal){
                p <- add_row_signal(p,
                                    signal[[i]], 
                                    name = signal_name)
              }
              if (summary){
                p <- suppressWarnings(add_col_summary(p,
                                    groups = TRUE, 
                                    showlegend = FALSE, 
                                    yname = summary_yaxis))
              }
              
              if (show_xlabels){
                p <- add_col_labels(p,
                                    ticktext = ticktext, 
                                    tickvals = tickvals)
              }
              if (!is.null(col_title)){
                p <- add_col_title(p, col_title)
              }
              if (!is.null(names(data)[i]))
                p <- add_col_title(p, names(data)[i], side = "top")
            }
            
            
            return(p)
            
          })          






default_end <- function(x){
  stopifnot(length(x) > 2)
  end = as.numeric(x[length(x)]) + as.numeric(x[2]) - as.numeric(x[1])
  as.character(end)
}



#' @export
#' @importFrom GenomicRanges distanceToNearest
add_dist_to_tss <- function(p,
                            ranges,
                            tss,
                            ...){
  
  d = log10(abs(mcols(distanceToNearest(ranges, tss, 
                                        ignore.strand = TRUE))$distance) + 1)
  
  
  p <- add_row_plot(p,
                    x = d,
                    side = "right", 
                    type = "scatter", 
                    mode = "markers",
                    x_layout = list(title = "Distance<br>to TSS<br>(log10)", 
                                    range = c(-0.5,6),
                                    tickvals = c(0,3,6),
                                    zeroline = FALSE), 
                    size = 0.3,
                    buffer = 0.04,
                    showlegend = FALSE)
  
  return(p)
}




