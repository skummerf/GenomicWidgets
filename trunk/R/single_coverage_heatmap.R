

# Plotting function.  Does not do any ordering. 

single_coverage_heatmap <- function(mat, 
                                    x = default_x(mat),
                                    y = default_y(mat),                    
                row_order = c("none","hclust","kmeans","groups","signal"),
                k = NULL,
                groups = NULL,
                clust_dist = stats::dist,
                signal = NULL,
                name = "Signal",
                summary = TRUE,
                source = "HM"){
  
  # TO DO: Add argument check
  
  row_order = match.arg(row_order)
  
  # 
  
  dendro = NULL
  
  if (row_order == "hclust"){
    dendro = flashClust::hclust(clust_dist(mat))
    row_order = dendro$order
    mat = mat[row_order,]
    if (!is.null(signal)){
      signal = signal[row_order]
    }
    if (!is.null(k)){
      groups = cutree(dendro, k = k)[row_order]
    }
    y = y[row_order]
  } else if (row_order == "kmeans"){
    stopifnot(!is.null(k))
    groups = kmeans(mat, centers = k)$cluster
    row_order = order(groups)
    groups = groups[row_order]
    mat = mat[row_order,]
    if (!is.null(signal)){
      signal = signal[row_order]
    }
    y = y[row_order]
  } else if (row_order == "groups"){
    row_order = order(groups)
    groups = groups[row_order]
    mat = mat[row_order,]
    if (!is.null(signal)){
      signal = signal[row_order]
    }
    y = y[row_order]
  } else if (row_order == "signal"){
    stopifnot(!is.null(signal))
    row_order = order(signal, decreasing = FALSE)
    groups = groups[row_order]
    mat = mat[row_order,]
    signal = signal[row_order]
    y = y[row_order]
  } else{
    row_order = 1:nrow(mat)
  }
  
  p <- function(){
    single_coverage_heatmap_helper(mat,
                                      x,
                                      y,
                                      groups,
                                      dendro,
                                      summary,
                                      signal,
                                      name,
                                      source)
  }
  return(list(plot = p, row_order = row_order, dendro = dendro))
}

default_x <- function(mat){
  if (is.null(colnames(mat))){
    return(1:ncol(mat))
  } else{
    colnames(mat)
  }
}

default_y <- function(mat){
  if (is.null(rownames(mat))){
    return(1:nrow(mat))
  } else{
    rownames(mat)
  }
}


single_coverage_heatmap_helper <- function(mat,
                                    x = default_x(mat),
                                    y = default_y(mat),
                                    groups = NULL,
                                    dendro = NULL,#hclust object
                                    summary = TRUE,
                                    signal = NULL,
                                    name = "Signal",
                                    source = "HM",
                                    zmin = min(mat),
                                    zmax = max(mat)){
  
  
  out<- plot_ly(z = mat, 
                x = x, 
                #text = matrix(y, nrow = nrow(mat), ncol = ncol(mat), byrow = FALSE),
                #hoverinfo="x+y+z",
                type="heatmap",
                colorbar = list(title = name,
                                     len = 0.3, 
                                     y = 0.6,
                                     x = 1.05), 
                xaxis = "x", 
                yaxis = "y",
                zmin = zmin,
                zmax = zmax,
                source = source) 

  # establish layout
  
  l <- list(dendro = !is.null(dendro),
            groups = !is.null(groups),
            summary = summary,
            signal = !is.null(signal))


  # groups? x3
  
  if (l$groups){
    if (!is.factor(groups)) groups <- as.factor(groups)
    ngroups = length(levels(groups))
    out <- out %>% add_trace(z = matrix(as.numeric(groups), ncol = 1), 
                             type = "heatmap",
                          colorscale = dcolorscale(ngroups), 
                          xaxis = "x3",
                          colorbar = list(x = 1.05, 
                                          y = 0.93, 
                                          len = 0.3, 
                                          title = "Row\nClusters",
                                          ticktext = levels(groups),
                                          tickvals = 1:length(levels(groups))))
  }
  
  # summary? y2
  
  if (l$summary){
    if (l$groups){
      out <- out %>% add_trace(x = x, #rep(0:(ncol(mat)-1),length(levels(groups))),
                               y = unlist(lapply(levels(groups),function(x) colMeans(mat[which(groups == x),]))),
                               color = rep(levels(groups), each = ncol(mat)),
                               yaxis = "y2",
                               mode="lines",
                               showlegend = FALSE)
    } else{
      out <- out %>% add_trace(x = x,#0:(ncol(mat)-1), 
                               y = colMeans(mat),
                               yaxis = "y2",
                               mode="lines",
                               showlegend = FALSE)
    }
  }
  
  # signal_strength? x4

 if (l$signal){
  out <- out %>% add_trace(z = matrix(signal, ncol = 1), 
                           type = "heatmap",
                           colorscale = "Reds", 
                           xaxis = "x4",
                           colorbar = list(x = 1.05, 
                                           y = 0.27, 
                                           len = 0.3, 
                                           title = "Total Signal"))
  
 }
  
  # dendro? x2
  
  if (l$dendro){
    #Construct dendrogram
    dendro_row = ggdendro::dendro_data(dendro)
    dendro_segments_row = dendro_row$segments
    dendro_segments_row[,c("y","yend")] = (dendro_segments_row[,c("y","yend")])*-1 
    dendro_df = rbind(data.frame(x = dendro_segments_row$y, y = dendro_segments_row$x - 1,group = 1:nrow(dendro_segments_row)),
                      data.frame(x = dendro_segments_row$yend, y = dendro_segments_row$xend - 1,group =1:nrow(dendro_segments_row)))
    
    out <- out %>% add_trace(x = dendro_df$x, 
                             y = dendro_df$y, 
                             group = dendro_df$group, 
                             mode = "lines", 
                             showlegend = F, 
                             line = list(color = "gray"),
                             xaxis = "x2",
                             hoverinfo = "none")
  }
  
  
  # put it together
  
  # x axis
  xm <- 0
  xm2 <- 1
  if (l$dendro){
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
  if (l$groups){
    out <- out %>% layout(xaxis3 = list(domain= c(ifelse(l$dendro,0.15,0),
                                                  ifelse(l$dendro,0.23,0.08)),
                                        anchor = "y",
                                        side = "bottom",
                                        tickvals = 0,
                                        ticktext = "Cluster",
                                        ticks = "",
                                        zeroline = FALSE,
                                        showline = FALSE,
                                        showgrid = FALSE))
    
    xm <- xm + 0.1
  } 
  if (l$signal){
    out <- out %>% layout(xaxis4 = list(domain = c(0.92,1),
                                        anchor = "y",
                                        side = "bottom",
                                        tickvals = 0,
                                        ticktext = "Signal",
                                        ticks = "",
                                        zeroline = FALSE,
                                        showline = FALSE,
                                        showgrid = FALSE))
    xm2 <- 0.9
  }
  
  #yaxis
  ym = 1
  if (l$summary){
    out <- out %>% layout(yaxis2 = list(domain = c(0.8,1),
                                        anchor = "x"))#,
                                        #zeroline = FALSE,
                                        #showline = FALSE))
    ym = 0.77
  }
  
  
  #main
  out <- out %>% layout(xaxis = list(domain = c(xm,xm2),
                                     anchor = "y",
                                     side = "bottom",
                                     tickangle = -90,
                                     #tickvals = 0:(ncol(mat)-1),
                                     #ticktext = x,
                                     ticks = "",
                                     zeroline = FALSE,
                                     showline = FALSE,
                                     showgrid = FALSE),
                        yaxis = list(domain = c(0,ym),
                                     anchor = "x",
                                     tickvals = 0:(nrow(mat)-1),
                                     ticktext = y,
                                     ticks = "",
                                     zeroline = FALSE,
                                     showline = FALSE,
                                     showgrid = FALSE,
                                     showticklabels = FALSE),
                        hovermode = 'closest')
    
  return(out)
  
}

  