default_end <- function(x){
  stopifnot(length(x) > 2)
  end = as.numeric(x[length(x)]) + as.numeric(x[2]) - as.numeric(x[1])
  as.character(end)
}

get_left_row_groups <- function(p){
  stopifnot(p$orientation == "horizontal")
  xcandidates_id <- p$plots %>% filter(type == "row_groups") %>% select(xid) %>% unlist()
  xcandidates <- p$xaxes %>% filter(xid %in% xcandidates_id)
  min_start <- min(xcandidates$domain_start)
  xaxis <- xcandidates %>% filter(domain_start == min_start) %>% select(xid) %>% unlist()
  gr_pid <- p$plots %>% filter(xid == xaxis, type == "row_groups") %>% select(pid) %>% unlist()
  groups <- p$plots %>% filter(pid == gr_pid) %>% 
    select(input_data) %>% 
    unlist(recursive = FALSE)
  return(groups[[1]])
}

#' single_coverage_heatmap
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
#' @return iheatmap object
#' @export
#' @author Alicia Schep
single_coverage_heatmap <- function(mat, 
                                    x = iheatmapr:::default_x(mat),
                                    y = iheatmapr:::default_y(mat),
                                    include = 1000,
                                    include_method = c("signal","first","random"),
                                    row_order = c("none","hclust","kmeans","groups","signal"),
                                    k = NULL,
                                    groups = NULL,
                                    clust_dist = stats::dist,
                                    signal = rowMeans(mat, na.rm = TRUE),
                                    plot_signal = TRUE,
                                    name = "Coverage",
                                    signal_name = "Avg. Coverage",
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
  
  if (length(include) > 1){
    if (is.numeric(include)){
      keep <- include
    } else {
      keep <- which(y %in% include)
      if (length(keep) == 0) stop("Values for incluce don't match y")
    }
  } else if (include < nrow(mat)){
    if (include_method  == "signal"){
      keep <- which(signal >= quantile(signal, (length(signal) - include) / length(signal), na.rm = TRUE))
    } else if (include_method == "first"){
      keep <- 1:include
    } else if (include_method == "random"){
      keep <- sample(1:nrow(mat), include)
    }
  } else{
    keep <- 1:nrow(mat)
  }
  
  if (scale_method != "none") mat <- normalize_coverage_matrix(mat, method = scale_method, pct = pct, scalar = scale_factor)
  
  dendro = NULL
  
  if (row_order == "hclust"){
    mat_sub <- mat[keep,]
    dendro = fastcluster::hclust(clust_dist(mat_sub))
    row_order = keep[dendro$order]
    if (!is.null(k)){
      groups = rep(NA, nrow(mat))
      groups[keep] = cutree(dendro, k = k)
    }
  } else if (row_order == "kmeans"){
    mat_sub <- mat[keep,]
    stopifnot(!is.null(k))
    groups = rep(NA, nrow(mat))
    groups[keep] = kmeans(mat_sub, centers = k)$cluster
    row_order = keep[order(groups[keep])]
  } else if (row_order == "groups"){
    stopifnot(!is.null(groups))
    groups_sub <- groups[keep]
    row_order = keep[order(groups_sub)]
  } else if (row_order == "signal"){
    signal_sub <- signal[keep]
    row_order = keep[order(signal_sub, decreasing = FALSE)]
  } else{
    row_order = keep
  }
  
  
  p <- iheatmap(mat, 
                x = x, 
                y = y, 
                row_order = row_order, 
                name = name,
                x_categorical = FALSE,
                font = font,
                ...) 
  
  if (isTRUE(ticktext)){
    if ("0" %in% x || 0 %in% x){
      ticktext = c(start, "0",end)
      tickvals = as.numeric(c(x[1], 0, x[length(x)]))
    } else{
      tickvals = c(0, ncol(mats[[1]]) - 1)
      ticktext = c(start, end)
      tickvals = as.numeric(c(x[1], x[length(x)]))
    }
    
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
#' @return iheatmap object
#' @export
#' @author Alicia Schep
add_coverage_heatmap <- function(p,
                                 mat, 
                                 x = iheatmapr:::default_x(mat),                    
                                 groups = p$row_groups,
                                 signal = rowMeans(mat, na.rm = TRUE),
                                 plot_signal = TRUE,
                                 name = "Coverage",
                                 signal_name = "Avg. Coverage",
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
  
  force(signal)
  if (length(signal) != nrow(mat)){
    stop("Invalid signal input.  Must be vector of length nrow(mat) or TRUE/FALSE")
  }
  
  if (scale_method != "none") mat <- normalize_coverage_matrix(mat, method = scale_method, pct = pct, scalar = scale_factor)
  
  
  dendro = NULL
  
  p <- p %>% add_iheatmap(mat, 
                          x = x, 
                          name = name,
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
#' @return iheatmap object
#' @export
#' @import iheatmapr
#' @author Alicia Schep
multi_coverage_heatmap <- function(mats, 
                                   x = iheatmapr:::default_x(mats[[1]]),
                                   y = iheatmapr:::default_y(mats[[1]]), 
                                   include = 1000,
                                   include_method = c("signal","first","random"),
                                   row_order = c("none","hclust","kmeans","groups","signal"),
                                   k = NULL,
                                   groups = NULL,
                                   clust_dist = stats::dist,
                                   cluster_by = c("first","all"), 
                                   signal = lapply(mats, rowMeans, na.rm = TRUE),
                                   plot_signal = TRUE, 
                                   name = "Coverage",
                                   signal_name = "Avg. Coverage",
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
  
  if (length(include) > 1){
    if (is.numeric(include)){
      keep <- include
    } else {
      keep <- which(y %in% include)
      if (length(keep) == 0) stop("Values for incluce don't match y")
    }
  } else if (include < nrow(mats[[1]])){
    if (include_method  == "signal"){
      tmp_signal <- Reduce("+",lapply(signal, function(z) z / sum(z)))
      keep <- which(tmp_signal >= quantile(tmp_signal, (length(tmp_signal) - include) / length(tmp_signal)))
    } else if (include_method == "first"){
      keep <- 1:include
    } else if (include_method == "random"){
      keep <- sample(1:nrow(mats[[1]]), include)
    }
    
  } else{
    keep <- 1:nrow(mats[[1]])
  }
  
  if (scale_method != "none"){
    mats <- normalize_coverage_matrix(mats, method = scale_method, pct = pct, scalar = scale_factor)
  } 
  
  
  dendro = NULL
  
  if (row_order == "hclust"){
    mats_sub <- lapply(mats, function(z) z[keep,])
    if (cluster_by == "first"){
      dendro = fastcluster::hclust(clust_dist(mats_sub[[1]]))
    } else{
      dendro = fastcluster::hclust(clust_dist(do.call(cbind,mats_sub)))
    }
    row_order = keep[dendro$order]
    if (!is.null(k)){
      groups = rep(NA, nrow(mats[[1]]))
      groups[keep] = cutree(dendro, k = k)
    }
  } else if (row_order == "kmeans"){
    mats_sub <- lapply(mats, function(z) z[keep,])
    stopifnot(!is.null(k))
    if (cluster_by == "first"){
      groups = rep(NA, nrow(mats[[1]]))
      groups[keep] = kmeans(mats_sub[[1]], centers = k)$cluster 
    } else{
      groups = rep(NA, nrow(mats[[1]]))
      groups[keep] = kmeans(do.call(cbind,mats_sub), centers = k)$cluster     
    }
    row_order = keep[order(groups[keep])]
  } else if (row_order == "groups"){
    groups_sub <- groups[keep]
    row_order = keep[order(groups_sub)]
  } else if (row_order == "signal"){
    signal_sub <- lapply(signal, function(z) z[keep])
    if (cluster_by == "first"){
      row_order = keep[order(signal_sub[[1]], decreasing = FALSE)]
    } else{
      tmp_signal <- Reduce("+",lapply(signal_sub, function(z) z / sum(z)))
      row_order = keep[order(tmp_signal, decreasing = FALSE)]
    }
  } else{
    row_order = keep
  }
  
  zmax <- max(sapply(mats, function(x) max(x[row_order,], na.rm = TRUE)))
  zmin <- min(sapply(mats, function(x) min(x[row_order,], na.rm = TRUE)))
  
  signal_zmax <- max(sapply(signal, function(x) max(x[row_order], na.rm = TRUE)))
  signal_zmin <- min(sapply(signal, function(x) min(x[row_order], na.rm = TRUE)))
  
  
  p <- iheatmap(mats[[1]], 
                x = x, 
                y = y, 
                row_order = row_order, 
                name = name,
                zmin = zmin,
                zmax = zmax,
                font = font,
                x_categorical = FALSE,
                ...) 
  
  if (isTRUE(ticktext)){
    if ("0" %in% x || 0 %in% x){
      ticktext = c(start, "0",end)
      tickvals = as.numeric(c(x[1], 0, x[length(x)]))
    } else{
      tickvals = c(0, ncol(mats[[1]]) - 1)
      ticktext = c(start, end)
      tickvals = as.numeric(c(x[1], x[length(x)]))
    }
    
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
    summary_yaxis = paste0("summary", max(p$plots$pid) + 1)
    p <- p %>% add_col_summary(groups, showlegend = FALSE, y_layout = list(font = font), 
                               yname = summary_yaxis)
  }  
  if (plot_signal){
    p <- p %>% add_row_signal(signal[[1]], 
                              name = signal_name,
                              zmin = signal_zmin, 
                              zmax = signal_zmax,
                              x_layout = list(font = font))
  }
  
  p <- p %>% add_x_axis_title(names(mats)[1], side = "top", font = font)
  
  if (length(mats) > 1){
    for (i in 2:length(mats)){
      p <- p %>% add_iheatmap(mats[[i]], 
                              x = x, 
                              show_colorbar = FALSE,
                              zmax = zmax,
                              zmin = zmin,
                              x_categorical = FALSE,
                              ...) 
      
      if (plot_signal){
        p <- p %>% add_row_signal(signal[[i]], 
                                  name = signal_name,
                                  zmin = signal_zmin, 
                                  zmax = signal_zmax,
                                  show_colorbar = FALSE)
      }
      if (summary){
        p <- p %>% add_col_summary(groups, showlegend = FALSE, yname = summary_yaxis)
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
  
  return(p)
}




#' add_multi_coverage_heatmap
#' @param mat coverage matrix
#' @param x x axis labels
#' @param y y axis labels
#' @param groups pre-determined groups for rows
#' @param signal signal along row
#' @param name name of colorbar
#' @param summary make summary plot, boolean
#' @param source source name in plotly
#' @param scale scale rows? or none
#' @param scale_method method to use for scaling
#' @return iheatmap object
#' @export
#' @import iheatmapr
#' @author Alicia Schep
add_multi_coverage_heatmap <- function(p,
                                       mats, 
                                       x = iheatmapr:::default_x(mats[[1]]),
                                       groups = NULL,
                                       signal = lapply(mats, rowMeans, na.rm = TRUE),
                                       plot_signal = TRUE, 
                                       name = "Coverage",
                                       signal_name = "Avg. Coverage",
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
  
  scale_method <- match.arg(scale_method)
  
  if (length(unique(lapply(mats, nrow))) > 1) stop("All input matrices must be of same length")
  
  force(signal)
  if (length(signal[[1]]) != nrow(mats[[1]])){
    stop("Invalid signal input.  Must be vector of length nrow(mat) or TRUE/FALSE")
  }
  
  if (scale_method != "none"){
    mats <- normalize_coverage_matrix(mats, method = scale_method, pct = pct, scalar = scale_factor)
  } 
  
  ro <- p$yaxes %>% filter(yid == "y") %>% select(row_order) %>% unlist()
  
  zmax <- max(sapply(mats, function(x) max(x[ro,], na.rm = TRUE)))
  zmin <- min(sapply(mats, function(x) min(x[ro,], na.rm = TRUE)))
  
  signal_zmax <- max(sapply(signal, function(x) max(x[ro], na.rm = TRUE)))
  signal_zmin <- min(sapply(signal, function(x) min(x[ro], na.rm = TRUE)))
  
  
  p <- p %>% add_iheatmap(mats[[1]], 
                          x = x, 
                          name = name,
                          zmin = zmin,
                          zmax = zmax,
                          x_categorical = FALSE,
                          ...) 
  
  if (isTRUE(ticktext)){
    if ("0" %in% x || 0 %in% x){
      ticktext = c(start, "0",end)
      tickvals = as.numeric(c(x[1], 0, x[length(x)]))
    } else{
      tickvals = c(0, ncol(mats[[1]]) - 1)
      ticktext = c(start, end)
      tickvals = as.numeric(c(x[1], x[length(x)]))
    }
    
    p <- p %>% add_x_axis_labels(ticktext = ticktext, tickvals = tickvals, font = font)
    
  }  else{
    tickvals <- NULL
  }
  
  if (is.character(xlab) && nchar(xlab) > 0){
    p <- p %>% add_x_axis_title(xlab, font = font)
  } else{
    xlab <- NULL
  }
  
  if (summary){
    summary_yaxis = paste0("summary", max(p$plots$pid) + 1)
    p <- p %>% add_col_summary(groups, showlegend = FALSE, y_layout = list(font = font),
                               yname = summary_yaxis)
  }  
  if (plot_signal){
    p <- p %>% add_row_signal(signal[[1]], 
                              name = signal_name,
                              zmin = signal_zmin, 
                              zmax = signal_zmax,
                              x_layout = list(font = font))
  }
  
  p <- p %>% add_x_axis_title(names(mats)[1], side = "top", font = font)
  
  if (length(mats) > 1){
    for (i in 2:length(mats)){
      if (share_z){
        p <- p %>% add_iheatmap(mats[[i]], 
                                x = x, 
                                show_colorbar = FALSE,
                                zmax = zmax,
                                zmin = zmin,
                                ...) 
      } else{
        p <- p %>% add_iheatmap(mats[[i]], 
                                x = x, 
                                name = names(mats)[i],
                                ...) 
      }
      
      
      if (plot_signal){
        p <- p %>% add_row_signal(signal[[i]], 
                                  name = signal_name,
                                  zmin = signal_zmin, 
                                  zmax = signal_zmax,
                                  show_colorbar = FALSE)
      }
      if (summary){
        p <- p %>% add_col_summary(groups, showlegend = FALSE, yname = summary_yaxis)
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




