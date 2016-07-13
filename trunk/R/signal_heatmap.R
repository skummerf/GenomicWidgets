


cluster_heatmap <- function(mat, 
                            x = default_x(mat),
                            y = default_y(mat),                   
                            row_order = c("none","hclust","kmeans","groups"),
                            col_order = c("none","hclust","kmeans","groups"),
                            row_groups = NULL,
                            col_groups = NULL,
                            row_anno = NULL,
                            col_anno = NULL,
                            row_k = NULL,
                            col_k = NULL,
                            row_clust_dist = stats::dist,
                            col_clust_dist = stats::dist,
                            name = "Signal",
                            summary = TRUE,
                            source = "HM",
                            scale = TRUE){
  
  # TO DO: Add argument check
  
  row_order = match.arg(row_order)
  col_order = match.arg(col_order)
  # 
  
  row_dendro = NULL
  col_dendro = NULL
  
  
  if (scale == TRUE){
    mat = t(scale(t(log10(mat+1))))
  }
  
  #Filter out rows with NA
  notna <- apply(mat,1, function(x) !any(is.na(x)))
  mat = mat[notna,]
  y = y[notna]
  if (!is.null(row_anno)){
    for ( i in seq_along(row_anno)){
      row_anno[[i]]$data = row_anno[[i]]$data[not_na]
    }
  }
  
  if (row_order == "hclust"){
    row_dendro = flashClust::hclust(row_clust_dist(mat))
    row_order = row_dendro$order
    mat = mat[row_order,]
    if (!is.null(row_anno)){
      for ( i in seq_along(row_anno)){
       row_anno[[i]]$data = row_anno[[i]]$data[row_order,]
      }
    }
    if (!is.null(row_k)){
      row_groups = cutree(row_dendro, k = row_k)[row_order]
    }
    y = y[row_order]
  } else if (row_order == "kmeans"){
    stopifnot(!is.null(row_k))
    row_groups = kmeans(mat, centers = row_k)$cluster
    row_order = order(row_groups)
    row_groups = row_groups[row_order]
    mat = mat[row_order,]
    if (!is.null(row_anno)){
      for ( i in seq_along(row_anno)){
        row_anno[[i]]$data = row_anno[[i]]$data[row_order,]
      }
    }
    y = y[row_order]
  } else if (row_order == "groups"){
    row_order = order(groups)
    row_groups = row_groups[row_order]
    mat = mat[row_order,]
    if (!is.null(row_anno)){
      for ( i in seq_along(row_anno)){
        row_anno[[i]]$data = row_anno[[i]]$data[row_order,]
      }
    }
    y = y[row_order]
  } else{
    row_order = 1:nrow(mat)
  }
  
  if (col_order == "hclust"){
    col_dendro = flashClust::hclust(col_clust_dist(t(mat)))
    col_order = col_dendro$order
    mat = mat[,col_order]
    if (!is.null(col_anno)){
      for ( i in seq_along(col_anno)){
        col_anno[[i]]$data = col_anno[[i]]$data[,col_order]
      }
    }
    if (!is.null(col_k)){
      col_groups = cutree(col_dendro, k = col_k)[col_order]
    }
    x = x[col_order]
  } else if (col_order == "kmeans"){
    stopifnot(!is.null(col_k))
    col_groups = kmeans(t(mat), centers = col_k)$cluster
    col_order = order(col_groups)
    col_groups = col_groups[col_order]
    mat = mat[,col_order]
    if (!is.null(col_anno)){
      for ( i in seq_along(col_anno)){
        col_anno[[i]]$data = col_anno[[i]]$data[,col_order]
      }
    }
    x = x[col_order]
  } else if (col_order == "groups"){
    col_order = order(groups)
    col_groups = col_groups[col_order]
    mat = mat[,col_order]
    if (!is.null(col_anno)){
      for ( i in seq_along(col_anno)){
        col_anno[[i]]$data = col_anno[[i]]$data[,col_order]
      }
    }
    x = x[col_order]
  } else{
    col_order = 1:ncol(mat)
  }  
  

  p <- function(){
    signal_heatmap_helper(mat = mat,
                          x = x,
                          y = y,
                          row_groups = row_groups,
                          col_groups = col_groups,
                          row_dendro = row_dendro,
                          col_dendro = col_dendro,
                          row_anno = row_anno,
                          col_anno = col_anno,
                          name = name,
                          source = source)
  }
  return(list(plot = p, 
              row_order = row_order, 
              row_dendro = row_dendro,
              col_order = col_order,
              col_dendro = col_dendro))
  
}


signal_heatmap_helper <- function(mat,
                                  x = default_x(mat),
                                  y = default_y(mat),
                                  row_groups = NULL,
                                  col_groups = NULL,
                                  row_dendro = NULL,#hclust object
                                           col_dendro = NULL,
                                           row_anno = NULL,
                                           col_anno = NULL,
                                           name = "Signal",
                                           source = "HM",
                                           zmin = min(mat),
                                           zmax = max(mat)){
  
  
  out<- plot_ly(z = mat, 
                type="heatmap",
                colorbar = list(title = name,
                                len = 0.3, 
                                y = 0.93,
                                thickness = 15,
                                x = 1.05), 
                xaxis = "x", 
                yaxis = "y",
                zmin = zmin,
                zmax = zmax,
                source = source) 
  
  # establish layout
  
  l <- list(row_dendro = !is.null(row_dendro),
            row_groups = !is.null(row_groups),
            row_anno = !is.null(row_anno),
            col_dendro = !is.null(col_dendro),
            col_groups = !is.null(col_groups),
            col_anno = !is.null(col_anno))
  
  
  # row dendro? x2
  
  if (l$row_dendro){
    #Construct dendrogram
    dendro_row = ggdendro::dendro_data(row_dendro)
    dendro_segments_row = dendro_row$segments
    dendro_segments_row[,c("y","yend")] = (dendro_segments_row[,c("y","yend")])*-1 
    row_dendro_df = rbind(data.frame(x = dendro_segments_row$y, y = dendro_segments_row$x - 1,group = 1:nrow(dendro_segments_row)),
                      data.frame(x = dendro_segments_row$yend, y = dendro_segments_row$xend - 1,group =1:nrow(dendro_segments_row)))
    
    out <- out %>% add_trace(x = row_dendro_df$x, 
                             y = row_dendro_df$y, 
                             group = row_dendro_df$group, 
                             mode = "lines", 
                             showlegend = F, 
                             line = list(color = "gray"),
                             xaxis = "x2",
                             hoverinfo = "none")
  }
  
  # col dendro? y2
  
  if (l$col_dendro){
    #Construct dendrogram
    dendro_col = ggdendro::dendro_data(col_dendro)
    dendro_segments_col = dendro_col$segments
    col_dendro_df = rbind(data.frame(y = dendro_segments_col$y, x = dendro_segments_col$x - 1,group = 1:nrow(dendro_segments_col)),
                      data.frame(y = dendro_segments_col$yend, x = dendro_segments_col$xend - 1,group =1:nrow(dendro_segments_col)))
    
    out <- out %>% add_trace(x = col_dendro_df$x, 
                             y = col_dendro_df$y, 
                             group = col_dendro_df$group, 
                             mode = "lines", 
                             showlegend = F, 
                             line = list(color = "gray"),
                             yaxis = "y2",
                             hoverinfo = "none")
  }
  
  # row groups? x3
  
  if (l$row_groups){
    if (!is.factor(row_groups)) row_groups <- as.factor(row_groups)
    nrowgroups = length(levels(row_groups))
    out <- out %>% add_trace(z = matrix(as.numeric(row_groups), ncol = 1), 
                             type = "heatmap",
                             colorscale = dcolorscale(nrowgroups), 
                             xaxis = "x3",
                             colorbar = list(x = 1.05, 
                                             y = 0.6, 
                                             len = 0.3, 
                                             thickness = 15,
                                             title = "Rows",
                                             ticktext = levels(row_groups),
                                             tickvals = 1:nrowgroups))
  }
  
  # col groups? y3
  
  if (l$col_groups){
    if (!is.factor(col_groups)) col_groups <- as.factor(col_groups)
    ncolgroups = length(levels(col_groups))
    out <- out %>% add_trace(z = matrix(as.numeric(col_groups), nrow = 1), 
                             type = "heatmap",
                             colorscale = dcolorscale(ncolgroups), 
                             yaxis = "y3",
                             colorbar = list(x = 1.05, 
                                             y = 0.27, 
                                             len = 0.3, 
                                             thickness = 15,
                                             title = "Cols",
                                             ticktext = levels(col_groups),
                                             tickvals = 1:ncolgroups))
  }
  
  
  # Row Anno? x4 +
  
  if (l$row_anno){
    stopifnot(is.list(row_anno))
    for ( i in seq_along(row_anno)){
      cbx <- 1.25 + (i %/% 3)*0.15
      cby <- 0.93 - ((i %% 3)-1)*0.33
      if (!inherits(row_anno[[i]]$data,"matrix")){
        tmpz <- matrix(row_anno[[i]]$data, nrow = nrow(mat)) 
      } else{
        tmpz <- row_anno[[i]]$data
      }
      out <- do.call(add_trace, list(out, z = tmpz, 
                                     type = "heatmap",
                                     colorscale = row_anno[[i]]$colorscale, 
                                     xaxis = paste0("x",3 + i),
                                     yaxis = "y",
                                     colorbar = list(x = cbx, 
                                                     y = cby, 
                                                     len = 0.3, 
                                                     thickness = 15,
                                                     title = names(row_anno)[i])))
    }
  }
  
  # Col Anno? y4 +
  
  if (l$col_anno){
    stopifnot(is.list(col_anno))
    for ( i in seq_along(col_anno)){
      cbx <- 1.25 + (length(col_anno) + (i %/% 3))*0.15
      cby <- 0.93 - ((i %% 3)-1)*0.33
      if (!inherits(col_anno[[i]]$data,"matrix")){
        tmpz <- matrix(col_anno[[i]]$data, ncol = ncol(mat)) 
      } else{
        tmpz <- col_anno[[i]]$data
      }
      out <- do.call(add_trace, list(out, 
                                     z = tmpz, 
                               type = "heatmap",
                               colorscale = col_anno[[i]]$colorscale, 
                               yaxis = paste0("y",3 + i),
                               xaxis = "x",
                               colorbar = list(x = cbx, 
                                               y = cby, 
                                               len = 0.3, 
                                               thickness = 15,
                                               title = names(col_anno)[i])))
    }
  }
  
  
  # put it together
  
  # x axis
  xm <- 0
  xm2 <- 1
  if (l$row_dendro){
    out <- out %>% layout(xaxis2 = list(domain= c(0,0.13),
                                        anchor = "y",
                                        title = "",
                                        zeroline = FALSE,
                                        showline = FALSE,
                                        showticklabels = FALSE,
                                        showgrid = FALSE,
                                        ticks = ""))
    
    xm <- 0.15
  } 
  if (l$row_groups){
    out <- out %>% layout(xaxis3 = list(domain= c(ifelse(l$row_dendro,0.15,0),
                                                  ifelse(l$row_dendro,0.2,0.05)),
                                        anchor = "y",
                                        side = "bottom",
                                        showticklabels = FALSE,
                                        tickvals = 0,
                                        ticktext = "Row Clusters",
                                        ticks = "",
                                        zeroline = FALSE,
                                        showline = FALSE,
                                        showgrid = FALSE))
    
    xm <- xm + 0.07
  } 
  
  
  if (l$row_anno){
    
    row_anno_layout <- vector(mode="list", length=length(row_anno))
    names(row_anno_layout) <- sapply(1:length(row_anno), function(x) paste0("xaxis", x +3)) 
    
    for (i in 1:length(row_anno)){
      row_anno_layout[[i]] = list(domain = c(xm2 - 0.05,xm2),
                                  anchor = "y",
                                  side = "bottom",
                                  tickangle = -90,
                                  tickvals = 0,
                                  ticktext = names(row_anno)[i],
                                  ticks = "",
                                  zeroline = FALSE,
                                  showline = FALSE,
                                  showgrid = FALSE)
      xm2 <- xm2 - 0.07
    }
    out <- do.call(layout, c(list(out), row_anno_layout))
  }
  
  ym = 0
  ym2 = 1
  
  if (l$col_dendro){
    out <- out %>% layout(yaxis2 = list(domain= c(0.87,1),
                                        anchor = "x",
                                        title = "",
                                        zeroline = FALSE,
                                        showline = FALSE,
                                        showticklabels = FALSE,
                                        showgrid = FALSE,
                                        ticks = ""))
    
    ym2 <- 0.85
  } 
  if (l$col_groups){
    out <- out %>% layout(yaxis3 = list(domain= c(ifelse(l$col_dendro,0.8,0.95),
                                                  ifelse(l$col_dendro,0.85,1)),
                                        anchor = "x",
                                        side = "bottom",
                                        tickvals = 0,
                                        ticktext = "Column Clusters",
                                        ticks = "",
                                        showticklabels = FALSE,
                                        zeroline = FALSE,
                                        showline = FALSE,
                                        showgrid = FALSE))
    
    ym2 <- ym2 - 0.07
  } 
  
  
  if (l$col_anno){
    
    col_anno_layout <- vector(mode="list", length=length(col_anno))
    names(col_anno_layout) <- sapply(1:length(col_anno), function(x) paste0("yaxis", x +3)) 
    
    for (i in seq_along(col_anno)){
      col_anno_layout[[i]] = list(domain = c((i-1)*0.07, (i-1)*0.07+ 0.05),
                                  anchor = "x",
                                  #side = "bottom",
                                  tickvals = 0,
                                  ticktext = names(col_anno)[i],
                                  ticks = "",
                                  zeroline = FALSE,
                                  showline = FALSE,
                                  showgrid = FALSE)
      
    }
    out <- do.call(layout, c(list(out), col_anno_layout))
    ym <- ym + 0.07*length(col_anno) + 0.05
  }
  
  

  
  #main
  out <- out %>% layout(xaxis = list(domain = c(xm,xm2),
                                     anchor = ifelse(l$col_anno, paste0("y", 3 + length(col_anno)),"x"),
                                     side = "bottom",
                                     title = "",
                                     tickvals = 0:(ncol(mat)-1),
                                     ticktext = x,
                                     ticks = "",
                                     zeroline = FALSE,
                                     showline = FALSE,
                                     showgrid = FALSE),
                        yaxis = list(domain = c(ym,ym2),
                                     anchor = ifelse(l$row_anno, paste0("x", 3 + length(row_anno)),"x"),
                                     tickvals = 0:(nrow(mat)-1),
                                     ticktext = y,
                                     ticks = "",
                                     tickangle = -90,
                                     zeroline = FALSE,
                                     showline = FALSE,
                                     showgrid = FALSE,
                                     showticklabels = FALSE))
  
  return(out)
  
}
