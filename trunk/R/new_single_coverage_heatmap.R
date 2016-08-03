default_end <- function(x){
  stopifnot(length(x) > 2)
  end = as.numeric(x[length(x)]) + as.numeric(x[2]) - as.numeric(x[1])
  as.character(end)
}

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
#' @return iHeatmap object
#' @export
#' @author Alicia Schep
new_single_coverage_heatmap <- function(mat, 
                                    x = iHeatmapR:::default_x(mat),
                                    y = iHeatmapR:::default_y(mat),
                                    include = 1000,
                                    include_method = c("signal","first","random"),
                                    row_order = c("none","hclust","kmeans","groups","signal"),
                                    k = NULL,
                                    groups = NULL,
                                    clust_dist = stats::dist,
                                    signal = rowSums(mat, na.rm = TRUE),
                                    plot_signal = TRUE,
                                    name = "Coverage",
                                    signal_name = "Aggregate",
                                    summary = TRUE,
                                    source = "HM",
                                    scale_method = c("localRms", 
                                                     "localMean", 
                                                     "localNonZeroMean", 
                                                     "PercentileMax", 
                                                     "scalar", 
                                                     "none"),
                                    pct = 0.95,
                                    scale_factor = 1, 
                                    ticktext = TRUE,
                                    start = x[1],
                                    end = default_end(x),     
                                    xlab = "Position",
                                    font = list(size = 8),
                                    ...){
  
  # TO DO: Add argument check
  
  row_order = match.arg(row_order)
  scale_method <- match.arg(scale_method)
  include_method <- match.arg(include_method)
  
  force(signal)
  if (length(signal) != nrow(mat)){
    stop("Invalid signal input.  Must be vector of length nrow(mat) or TRUE/FALSE")
  }
  
  if (include < nrow(mat)){
    if (include_method  == "signal"){
      keep <- which(signal >= quantile(signal, (length(signal) - include) / length(signal), na.rm = TRUE))
    } else if (include_method == "first"){
      keep <- 1:include
    } else if (include_method == "random"){
      keep <- sample(1:nrow(mat), include)
    }
    mat <- mat[keep,]
    y <- y[keep]
    signal <- signal[keep]
  }
  
  if (scale_method != "none") mat <- normalize_coverage_matrix(mat, method = scale_method, pct = pct, scalar = scale_factor)
  
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
    row_order = order(signal, decreasing = FALSE)
  } else{
    row_order = 1:nrow(mat)
  }
  
  p <- make_main_hm(mat, 
                    x = x, 
                    y = y, 
                    row_order = row_order, 
                    name = name,
                    ...) 
  
  if (isTRUE(ticktext)){
    if ("0" %in% x){
      tickvals = c(0, which(x == "0") - 1, ncol(mat) - 1)
      ticktext = c(start, "0",end)
    } else{
      tickvals = c(0, ncol(mat) - 1)
      ticktext = c(start, end)
    }
  
    p <- p %>% add_x_axis_labels(ticktext = ticktext, tickvals = tickvals, font = font)
  
  } else if (!is.null(ticktext) && ticktext != FALSE){
    tickvals = which(x %in% ticktext) - 1
    p <- p %>% add_x_axis_labels(ticktext = ticktext, tickvals = tickvals, font = font)
  }
  
  if (is.character(xlab) && nchar(xlab) > 0){
    p <- p %>% add_x_axis_title(xlab, font = font)
  }
  
  if (!is.null(groups)){
    p <- p %>% add_row_groups(groups, side = "left")
  }
  if (!is.null(dendro)){
    p <- p %>% add_row_dendro(dendro, side = "left")
  }    
  if (plot_signal){
    p <- p %>% add_row_signal(signal, signal_name, 
                              x_layout = list(font = font))
  }
  if (summary){
    p <- p %>% add_col_summary(groups = groups, showlegend = FALSE,
                               y_layout = list(font = font))
  }
  
  p$row_groups <- groups
  
  p$layout$font <- font
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
                                 include = 1000,
                                 include_method = c("signal","first","random"),
                                 signal = rowSums(mat),
                                 plot_signal = TRUE,
                                 name = "Coverage",
                                 signal_name = "Aggregate",
                                 summary = TRUE,
                                 scale_method = c("localRms", 
                                                  "localMean", 
                                                  "localNonZeroMean", 
                                                  "PercentileMax", 
                                                  "scalar", 
                                                  "none"),
                                 pct = 0.95,
                                 scale_factor = 1, 
                                 ticktext = TRUE,
                                 start = x[1],
                                 end = default_end(x),     
                                 xlab = "Position",
                                 font = list(size = 8),
                                 ...){
  
  
  scale_method <- match.arg(scale_method)
  include_method <- match.arg(include_method)
  
  
  force(signal)
  if (length(signal) != nrow(mat)){
    stop("Invalid signal input.  Must be vector of length nrow(mat) or TRUE/FALSE")
  }
  
  if (include < nrow(mat)){
    if (include_method  == "signal"){
      keep <- which(signal >= quantile(signal, (length(signal) - include) / length(signal)))
    } else if (include_method == "first"){
      keep <- 1:include
    } else if (include_method == "random"){
      keep <- sample(1:nrow(mat), include)
    }
    mat <- mat[keep,]
    y <- y[keep]
    signal <- signal[keep]
  }
  

  if (scale_method != "none") mat <- normalize_coverage_matrix(mat, method = scale_method, pct = pct, scalar = scale_factor)
  
  
  dendro = NULL
  
  p <- p %>% add_main_hm(mat, 
                    x = x, 
                    name = name,
                    y_pos = "none",
                    x_pos = "none",
                    buffer = 0.1,
                    ...) 
  
  if (isTRUE(ticktext)){
    if ("0" %in% x){
      tickvals = c(0, which(x == "0") - 1, ncol(mat) - 1)
      ticktext = c(start, "0",end)
    } else{
      tickvals = c(0, ncol(mat) - 1)
      ticktext = c(start, end)
    }
    
    p <- p %>% add_x_axis_labels(ticktext = ticktext, tickvals = tickvals,
                                 font = font)
    
  } else if (!is.null(ticktext) && ticktext != FALSE){
    tickvals = which(x %in% ticktext) - 1
    p <- p %>% add_x_axis_labels(ticktext = ticktext, 
                                 tickvals = tickvals,
                                 font = font)
  }
  
  if (is.character(xlab) && nchar(xlab) > 0){
    p <- p %>% add_x_axis_title(xlab)
  }
  
  
  if (plot_signal){
    p <- p %>% add_row_signal(signal, signal_name, 
                              x_layout = list(font = font))
  }
  if (summary){
    p <- p %>% add_col_summary(groups, showlegend = FALSE, y_layout = list(font = font))
  }
  
  p$layout$font <- font
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
                                   include = 1000,
                                   include_method = c("signal","first","random"),
                                        row_order = c("none","hclust","kmeans","groups","signal"),
                                        k = NULL,
                                        groups = NULL,
                                        clust_dist = stats::dist,
                                        cluster_by = c("first","all"), 
                                        signal = lapply(mats, rowSums),
                                   plot_signal = TRUE, 
                                   name = "Coverage",
                                   signal_name = "Aggregate",
                                        summary = TRUE,
                                        source = "HM",
                                   scale_method = c("localRms", 
                                                    "localMean", 
                                                    "localNonZeroMean", 
                                                    "PercentileMax", 
                                                    "scalar", 
                                                    "none"),
                                   pct = 0.95,
                                   scale_factor = rep(1, length(mats)), 
                                   share_z = TRUE,
                                   ticktext = TRUE,
                                   start = x[1],
                                   end = default_end(x),     
                                   xlab = "Position",
                                   font = list(size = 8),
                                        ...){
  
  # TO DO: Add argument check
  
  row_order = match.arg(row_order)
  scale_method <- match.arg(scale_method)
  cluster_by <- match.arg(cluster_by)
  include_method <- match.arg(include_method)
  
  if (length(unique(lapply(mats, nrow))) > 1) stop("All input matrices must be of same length")
  
  force(signal)
  if (length(signal[[1]]) != nrow(mats[[1]])){
    stop("Invalid signal input.  Must be vector of length nrow(mat) or TRUE/FALSE")
  }
  
  if (include < nrow(mats[[1]])){
    if (include_method  == "signal"){
      tmp_signal <- Reduce("+",lapply(signal, function(z) z / sum(z)))
      keep <- which(tmp_signal >= quantile(tmp_signal, (length(tmp_signal) - include) / length(tmp_signal)))
    } else if (include_method == "first"){
      keep <- 1:include
    } else if (include_method == "random"){
      keep <- sample(1:nrow(mats[[1]]), include)
    }
    mats <- lapply(mats, function(z) z[keep,])
    y <- y[keep]
    signal <- lapply(signal, function(z) z[keep])
  }
  
  if (scale_method != "none"){
    mats <- normalize_coverage_matrix(mats, method = scale_method, pct = pct, scalar = scale_factor)
  } 
  
  if (share_z){
    zmax <- max(sapply(mats, max, na.rm = TRUE))
    zmin <- min(sapply(mats, min, na.rm = TRUE))
  }
  
  dendro = NULL
  
  if (row_order == "hclust"){
    if (cluster_by == "first"){
      dendro = flashClust::hclust(clust_dist(mats[[1]]))
    } else{
      dendro = flashClust::hclust(clust_dist(do.call(cbind,mats)))
    }
    row_order = dendro$order
    if (!is.null(k)){
      groups = cutree(dendro, k = k)
    }
  } else if (row_order == "kmeans"){
    stopifnot(!is.null(k))
    if (cluster_by == "first"){
      groups = kmeans(mats[[1]], centers = k)$cluster 
    } else{
       groups = kmeans(do.call(cbind,mats), centers = k)$cluster     
    }
    row_order = order(groups)
  } else if (row_order == "groups"){
    row_order = order(groups)
  } else if (row_order == "signal"){
    if (cluster_by == "first"){
      row_order = order(signal[[1]], decreasing = FALSE)
    } else{
      tmp_signal <- Reduce("+",lapply(signal, function(z) z / sum(z)))
      row_order = order(tmp_signal, decreasing = FALSE)
    }
  } else{
    row_order = 1:nrow(mats[[1]])
  }
  
  if (share_z){
    p <- make_main_hm(mats[[1]], 
                      x = x, 
                      y = y, 
                      row_order = row_order, 
                      name = name,
                      zmin = zmin,
                      zmax = zmax,
                      ...) 
  } else {
    p <- make_main_hm(mats[[1]], 
                    x = x, 
                    y = y, 
                    row_order = row_order, 
                    name = names(mats)[1],
                    ...)     
  }
  
  if (isTRUE(ticktext)){
    if ("0" %in% x){
      tickvals = c(0, which(x == "0") - 1, ncol(mats[[1]]) - 1)
      ticktext = c(start, "0",end)
    } else{
      tickvals = c(0, ncol(mats[[1]]) - 1)
      ticktext = c(start, end)
    }
    
    p <- p %>% add_x_axis_labels(ticktext = ticktext, tickvals = tickvals, font = font)
    
  } else if (!is.null(ticktext) && ticktext != FALSE){
    tickvals = which(x %in% ticktext) - 1
    p <- p %>% add_x_axis_labels(ticktext = ticktext, tickvals = tickvals, font = font)
  } else{
    tickvals <- NULL
  }
  
  if (is.character(xlab) && nchar(xlab) > 0){
    p <- p %>% add_x_axis_title(xlab, font = font)
  } else{
    xlab <- NULL
  }

  if (summary){
    p <- p %>% add_col_summary(groups, showlegend = FALSE, y_layout = list(font = font))
  }  
  if (plot_signal){
    p <- p %>% add_row_signal(signal[[1]], paste(names(mats)[1],signal_name,sep = "<br>"), 
                              x_layout = list(font = font))
  }
  
  p <- p %>% add_x_axis_title(names(mats)[1], side = "top", font = font)
  
  if (length(mats) > 1){
    for (i in 2:length(mats)){
      if (share_z){
        p <- p %>% add_main_hm(mats[[i]], 
                               x = x, 
                               show_colorbar = FALSE,
                               zmax = zmax,
                               zmin = zmin,
                               ...) 
      } else{
        p <- p %>% add_main_hm(mats[[i]], 
                             x = x, 
                             name = names(mats)[i],
                             ...) 
      }
      
      
      if (plot_signal){
        p <- p %>% add_row_signal(signal[[i]], paste(names(mats)[i],signal_name,sep = "<br>"), 
                                  x_layout = list(font = font))
      }
      if (summary){
        summary_yaxis = gsub("yaxis","y",names(p$data$xaxis)[which(sapply(p$data$xaxis, 
                                                                          function(x) x$type)=="col_summary")])
        p <- p %>% add_col_summary(groups, showlegend = FALSE, yaxis = summary_yaxis)
      }
      
      if (!is.null(tickvals)){
        p <- p %>% add_x_axis_labels(ticktext = ticktext, tickvals = tickvals, font = font)
      }
      if (!is.null(xlab)){
        p <- p %>% add_x_axis_title(xlab, font = font)
      }
      p <- p %>% add_x_axis_title(names(mats)[i], side = "top", font = font)
    }
  }
  
  if (!is.null(groups)){
    p <- p %>% add_row_groups(groups, side = "left")
  }
  if (!is.null(dendro)){
    p <- p %>% add_row_dendro(dendro, side = "left")
  }    
  p$layout$font <- font
  return(p)
}


#' @export
add_dist_to_tss <- function(p,
                            ranges,
                            tss,
                            ...){
  
  d = log10(abs(mcols(GenomicRanges::distanceToNearest(ranges, tss, ignore.strand = TRUE))$distance) + 1)
  
  
  p <- p %>% add_row_plot(x = d,
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




