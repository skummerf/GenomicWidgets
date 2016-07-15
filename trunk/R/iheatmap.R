# Main Functions -----------------------------------------------------------------------

#' plot_iHeatmap
#' 
#' Function to plot iHeatmap objects
#' 
#' @param ihm
#' @seealso 
#' @export
plot_iHeatmap <- function(ihm){
  stopifnot(inherits(ihm, "iHeatmap"))
  do.call(layout, c(ihm$plot, ihm$layout))
}

#' make_main_hm
#' 
#' Plots initial heatmap, creates iHeatmap object
#' @param mat matrix
#' @param x x axis labels (by default rownames of mat)
#' @param y y axis labels (by default colnames of mat)
#' @param row_order row order-- will be used for all subsequent elements sharing y axis
#' @param col_order column ordering for this heatmap
#' @param cbx x position of colorbar
#' @param cby y position of colorbar
#' @param colorscale colorscale to use
#' @param cbl length of colorbar
#' @param name name of colorbar
#' @param source source name -- used for shiny
#' @param zmin minimum for colorscale
#' @param zmax maximum for colorscale
#' @seealso \code{\link{add_main_hm}} 
#' @export
make_main_hm <- function(mat, 
                         x = default_x(mat),
                         y = default_y(mat),
                         row_order = 1:nrow(mat), 
                         col_order = 1:ncol(mat), 
                         cbx = 1.04, 
                         cby = 0.94, 
                         colorscale = "Viridis", 
                         cbl = 0.3,
                         name = "Signal",
                         source = "HM",
                         zmin = min(mat),
                         zmax = max(mat)){
  
  ## TO DO:  ADD ARGUMENT CHECKS!
  
  out <- list()
  class(out) <- "iHeatmap"
  
  out$row_order = row_order
  
  out$plot <- plot_ly(z = mat[row_order, col_order], 
                      type="heatmap",
                      colorscale = colorscale, 
                      colorbar = list(title = name,
                                      len = cbl, 
                                      y = cby,
                                      x = cbx), 
                      xaxis = "x", 
                      yaxis = "y",
                      source = source) 
  
  out$layout <- list(xaxis = list(domain = c(0,1),
                                  anchor = "y",
                                  side = "bottom",
                                  tickangle = -90,
                                  tickvals = 0:(ncol(mat)-1),
                                  ticktext = x,
                                  ticks = "",
                                  zeroline = FALSE,
                                  showline = FALSE,
                                  showgrid = FALSE),
                     yaxis = list(domain = c(0,1),
                                  anchor = "x",
                                  tickvals = 0:(nrow(mat)-1),
                                  ticktext = y,
                                  ticks = "",
                                  zeroline = FALSE,
                                  showline = FALSE,
                                  showgrid = FALSE,
                                  showticklabels = FALSE),
                     hovermode = 'closest')
  
  out$nx = 1
  out$ny = 1
  out$ncb = 1
  #out$curve_map = list()
  out$current_x = "x"
  out$xy = list(xaxis = c("yaxis"))
  
  return(out)
  
}


#' add_main_hm
#' 
#' Adds a "main heatmap" to an iHeatmap object
#' @param mat matrix
#' @param x x axis labels (by default rownames of mat)
#' @param y y axis labels (by default colnames of mat)
#' @param col_order column ordering for this heatmap
#' @param cbx x position of colorbar
#' @param cby y position of colorbar
#' @param colorscale colorscale to use
#' @param cbl length of colorbar
#' @param name name of colorbar
#' @param source source name -- used for shiny
#' @param zmin minimum for colorscale
#' @param zmax maximum for colorscale
#' @seealso \code{\link{make_main_hm}} 
#' @export
add_main_hm <- function(p,
                        mat, 
                        x = default_x(mat), 
                        y = default_y(mat),
                        row_order = 1:nrow(mat), 
                        col_order = 1:ncol(mat), 
                        xaxis = paste0("x",p$nx + 1), 
                        cbx = 1.04 + ((p$ncb-1) %/% 3) * 0.16, 
                        cby = 0.94 - ((p$ncb-1) %% 3)*0.32, 
                        colorscale = "Viridis", 
                        cbl = 0.3,
                        name = "Signal",
                        source = "HM",
                        size = 1,
                        show_colorbar = TRUE,
                        zmin = min(mat),
                        zmax = max(mat)){
  
  stopifnot(inherits(p, "iHeatmap"))
  p$plot <- p$plot %>% add_trace(z = mat[p$row_order, col_order], 
                                 type="heatmap",
                                 colorscale = colorscale, 
                                 zmin = zmin,
                                 zmax = zmax,
                                 colorbar = list(title = name,
                                                 len = cbl, 
                                                 y = cby,
                                                 x = cbx), 
                                 showscale = show_colorbar,
                                 xaxis = xaxis, 
                                 yaxis = "y",
                                 source = source) 
  
  p$layout = adjust_x_domains(p$layout, xaxis, size, p$nx, "right",
                              list(anchor = "y",
                                   side = "bottom",
                                   tickangle = -90,
                                   tickvals = 0:(ncol(mat)-1),
                                   ticktext = x,
                                   ticks = "",
                                   zeroline = FALSE,
                                   showline = FALSE,
                                   showgrid = FALSE))
  p$nx = p$nx + 1
  if (show_colorbar) p$ncb = p$ncb +1
  p$current_x = xaxis
  p$xy[[gsub("x","xaxis",xaxis)]] = c("yaxis")
  return(p)
}

#' add_row_dendro
#' 
#' Adds row dendrogram to iHeatmap object
#' @param mat matrix
#' @param dendro hclust object
#' @param xaxis xaxis to use
#' @param side side of plot on which to add dendro
#' @param size relative size of dendrogram (relative to the main heatmap)
#' @seealso \code{\link{make_main_hm}}, \code{\link{plot_iHeatmap}} 
#' @export
add_row_dendro <- function(p, 
                           dendro, 
                           xaxis = paste0("x",p$nx + 1),
                           side = c("left","right"),
                           size = 0.15){
  
  stopifnot(inherits(p, "iHeatmap"))
  side = match.arg(side)
  dendro_df <- dendro_to_segments(dendro, orientation = side)
  p$plot <- p$plot %>% add_trace(x = dendro_df$x, 
                                 y = dendro_df$y, 
                                 group = dendro_df$group, 
                                 mode = "lines", 
                                 showlegend = F, 
                                 line = list(color = "gray"),
                                 xaxis = xaxis,
                                 hoverinfo = "none")
  p$layout = adjust_x_domains(p$layout, xaxis, size, p$nx, side,
                              no_axis)
  p$nx = p$nx + 1
  
  return(p)
}

#' add_row_groups
#' 
#' Adds row groups to iHeatmap object
#' @param mat matrix
#' @param groups vector of group names
#' @param cbx x position of colorbar
#' @param cby y position of colorbar
#' @param colorscale colorscale to use
#' @param cbl length of colorbar
#' @param name name of colorbar
#' @param xaxis xaxis to use
#' @param side side of plot on which to add dendro
#' @param size relative size of dendrogram (relative to the main heatmap)
#' @param show_colorbar show the colorbar?
#' @seealso \code{\link{make_main_hm}}, \code{\link{plot_iHeatmap}}, \code{\link{add_col_groups}}
#' @export
add_row_groups <- function(p, 
                           groups, 
                           xaxis = paste0("x",p$nx + 1),
                           cbx = 1.04 + ((p$ncb-1) %/% 3) * 0.16, 
                           cby = 0.94 - ((p$ncb-1) %% 3)*0.32, 
                           colorscale = dcolorscale(length(levels(as.factor(groups)))),
                           side = c("left","right"),
                           name = "Rows", 
                           cbl = 0.3, 
                           size = 0.05,
                           show_colorbar = TRUE){
  side = match.arg(side)
  p$plot <- p$plot %>% add_trace(z = matrix(as.numeric(as.factor(groups))[p$row_order], ncol = 1), 
                                 text = groups,
                                 type = "heatmap",
                                 colorscale = colorscale, 
                                 xaxis = xaxis,
                                 showscale = show_colorbar,
                                 colorbar = list(x = cbx, 
                                                 y = cby, 
                                                 len = cbl, 
                                                 title = name))
  p$layout = adjust_x_domains(p$layout, 
                              xaxis = xaxis, 
                              size = size, 
                              nx = p$nx, 
                              side = side, 
                              new_x_layout = no_axis)
  p$nx = p$nx + 1
  if (show_colorbar) p$ncb = p$ncb +1
  return(p)
}



#' add_row_dendro
#' 
#' Adds row dendrogram to iHeatmap object
#' @param mat matrix
#' @param dendro hclust object
#' @param yaxis yaxis to use
#' @param xaxis xaxis to use
#' @param side side of plot on which to add dendro
#' @param size relative size of dendrogram (relative to the main heatmap)
#' @seealso \code{\link{make_main_hm}}, \code{\link{plot_iHeatmap}},  \code{\link{add_col_dendro}}
#' @export
add_col_dendro <- function(p, 
                           dendro, 
                           yaxis = paste0("y", p$ny + 1),
                           xaxis = p$current_x,
                           side = c("top","bottom"),
                           size = 0.15){
  
  side = match.arg(side)
  dendro_df <- dendro_to_segments(dendro, orientation = side)
  p$plot <- p$plot %>% add_trace(x = dendro_df$x, 
                                 y = dendro_df$y, 
                                 group = dendro_df$group, 
                                 mode = "lines", 
                                 showlegend = F, 
                                 line = list(color = "gray"),
                                 xaxis = xaxis,
                                 yaxis = yaxis,
                                 hoverinfo = "none")
  p$layout = adjust_y_domains(p$layout, 
                              yaxis = yaxis,
                              xaxis = xaxis, 
                              size = size, 
                              xy = p$xy, 
                              side = side,
                              new_y_layout = no_axis)
  p$ny = p$ny + 1
  p$xy[[gsub("x","xaxis",xaxis)]] = c(p$xy[[gsub("x","xaxis",xaxis)]], gsub("y","yaxis",yaxis))
  
  return(p)
}



#' add_col_groups
#' 
#' Adds column groups to iHeatmap object
#' @param mat matrix
#' @param groups vector of group names
#' @param yaxis yaxis to use
#' @param xaxis xaxis to use
#' @param cbx x position of colorbar
#' @param cby y position of colorbar
#' @param colorscale colorscale to use
#' @param cbl length of colorbar
#' @param name name of colorbar
#' @param side side of plot on which to add groups
#' @param size relative size of dendrogram (relative to the main heatmap)
#' @param show_colorbar show the colorbar?
#' @seealso \code{\link{make_main_hm}}, \code{\link{plot_iHeatmap}}, \code{\link{add_row_groups}}
#' @export
add_col_groups <- function(p, 
                           groups, 
                           yaxis = paste0("y",p$ny + 1),
                           xaxis = p$current_x,
                           cbx = 1.04 + ((p$ncb-1) %/% 3) * 0.16, 
                           cby = 0.94 - ((p$ncb-1) %% 3)*0.32, 
                           colorscale = dcolorscale(length(levels(as.factor(groups)))),
                           side = c("top","bottom"),
                           name = "Rows", 
                           cbl = 0.3, 
                           size = 0.05,
                           show_colorbar = TRUE){
  p$plot <- p$plot %>% add_trace(z = matrix(as.numeric(as.factor(groups)), nrow = 1), 
                                 text = groups,
                                 type = "heatmap",
                                 colorscale = colorscale, 
                                 xaxis = xaxis,
                                 yaxis = yaxis,
                                 colorbar = list(x = cbx, 
                                                 y = cby, 
                                                 len = cbl, 
                                                 title = name))
  p$layout <- adjust_y_domains(p$layout, 
                               yaxis = yaxis,
                               xaxis = xaxis, 
                               size = size, 
                               xy = p$xy, 
                               side = side,
                               new_y_layout = no_axis)
  p$ny <- p$ny + 1
  p$xy[[gsub("x","xaxis",xaxis)]] = c(p$xy[[gsub("x","xaxis",xaxis)]], gsub("y","yaxis",yaxis))
  p$ncb = p$ncb + 1
  
  return(p)
}

##TO DO:  Make functions for non-discrete row, col annotation

##TO DO:  Make functions for row/col summaries

# make_col_annot <- function(p, 
#                            anno,
#                            yaxis, 
#                            cbx, 
#                            cby, 
#                            colorscale, 
#                            name, 
#                            xaxis = "x", 
#                            cbl = 0.3){
#   if (!is.matrix(anno)){
#     anno =  matrix(anno, ncol = 1)
#   }
#   out <- p %>% add_trace(z = anno, 
#                          type = "heatmap",
#                          colorscale = colorscale, 
#                          xaxis = xaxis,
#                          yaxis = yaxis, 
#                          colorbar = list(x = cbx, 
#                                          y = cby, 
#                                          len = cbl, 
#                                          title = name))
#   return(out)
# }
# 
# make_row_annot <- function(p, anno, xaxis, cbx, cby, colorscale, cbl = 0.3){
#   out <- p %>% add_trace(z = matrix(anno, ncol = 1), 
#                          type = "heatmap",
#                          colorscale = colorscale, 
#                          xaxis = xaxis,
#                          colorbar = list(x = cbx, 
#                                          y = cby, 
#                                          len = cbl, 
#                                          title = name))
#   return(out)
# }
# 
# make_col_summary <- function(p, s, yaxis, xaxis = "x"){
#   out <- p %>% add_trace(s,
#                          yaxis = yaxis,
#                          xaxis = xaxis,
#                          showlegend = FALSE)
#   return(out)
# }
# 
# make_row_summary <- function(p, s, xaxis){
#   out <- p %>% add_trace(s,
#                          xaxis = xaxis,
#                          showlegend = FALSE)
#   return(out)
# }



# Utility Functions -------------------------------------------------------------------

'%ni%' = Negate('%in%')

dendro_to_segments <- function(dendro, orientation = c("left","bottom","right","top")){
  
  orientation <- match.arg(orientation)
  d <- ggdendro::dendro_data(dendro)$segments
  
  if (orientation == "left"){
    d[,c("y","yend")] = (d[,c("y","yend")])*-1 
    data = rbind(data.frame(x = d$y, y = d$x - 1,group = 1:nrow(d)),
                      data.frame(x = d$yend, y = d$xend - 1,group =1:nrow(d)))
    
  } else if (orientation == "top"){
    data = rbind(data.frame(x = d$x - 1, y = d$y,group = 1:nrow(d)),
                 data.frame(x = d$xend - 1, y = d$yend,group =1:nrow(d)))
  } else if (orientation == "right"){
    data = rbind(data.frame(x = d$y, y = d$x - 1,group = 1:nrow(d)),
                 data.frame(x = d$yend, y = d$xend - 1,group =1:nrow(d)))
    
  } else if (orientation == "bottom"){
    d[,c("y","yend")] = (d[,c("y","yend")])*-1 
    data = rbind(data.frame(x = d$x - 1, y = d$y,group = 1:nrow(d)),
                 data.frame(x = d$xend - 1, y = d$yend,group =1:nrow(d)))
  } 
  
  return(data)
}

row_summary <- function(mat, groups = NULL){
  
  if (is.null(groups)){
    out = data.frame(x = 0:(nrow(mat)-1), 
                     y = rowMeans(mat))
  } else{
    out = data.frame(x = rep(0:(nrow(mat)-1),length(levels(groups))),
                     y = unlist(lapply(levels(groups),function(x) rowMeans(mat[,which(groups == x)]))),
                     color = rep(levels(groups), each = nrow(mat)))
  }
  return(out)
}

col_summary <- function(mat, groups = NULL){
  if (is.null(groups)){
    out = data.frame(x = 0:(ncol(mat)-1), 
                     y = colMeans(mat))
  } else{
    out = data.frame(x = rep(0:(ncol(mat)-1),length(levels(groups))),
                     y = unlist(lapply(levels(groups),function(x) colMeans(mat[which(groups == x),]))),
                     color = rep(levels(groups), each = ncol(mat)))
  }
  return(out)
}

# Function to adjust x axis domains when adding in new xaxis
adjust_x_domains <- function(current_layout, 
                             xaxis, 
                             size, 
                             nx,
                             side = c("left","right"),
                             new_x_layout = list()){
  side = match.arg(side)
  current_x <- grep("xaxis", names(current_layout))
  x_domains <- lapply(current_x, function(x) current_layout[[x]]$domain)
  main_x_domain <- current_layout$xaxis$domain
  main_x_size <- diff(main_x_domain)
  x_sizes <- sapply(x_domains, diff)
  new_x_sizes <- c(x_sizes, main_x_size * size)
  buffer_size <-  0.04 / sqrt(nx)
  new_x_sizes <- new_x_sizes / (sum(new_x_sizes) + buffer_size* nx)
  x_order <- order(sapply(x_domains,min))
  if (side == "right"){
    j = 0
    for (i in seq_along(x_order)){
      current_layout[[current_x[x_order[i]]]]$domain = c(j,j + new_x_sizes[x_order[i]])
      j = j + new_x_sizes[x_order[i]] + buffer_size
    }
    current_layout[[gsub("x","xaxis",xaxis)]] = c(new_x_layout, list(domain = c(j,j + new_x_sizes[length(new_x_sizes)])))
  } else{
    j = new_x_sizes[length(new_x_sizes)] + buffer_size
    for (i in seq_along(x_order)){
      current_layout[[current_x[x_order[i]]]]$domain = c(j,j + new_x_sizes[x_order[i]])
      j = j + new_x_sizes[x_order[i]] + buffer_size
    }
    current_layout[[gsub("x","xaxis",xaxis)]] = c(new_x_layout, list(domain = c(0,new_x_sizes[length(new_x_sizes)])))
  }
  return(current_layout)
}

# Function to adjust y axis domains when adding in new yaxis
adjust_y_domains <- function(current_layout, 
                             yaxis,
                             xaxis,
                             size, 
                             xy,
                             side = c("top","bottom"),
                             new_y_layout = list()){
  side = match.arg(side)
  current_y <-  which(names(current_layout) %in% xy[[gsub("x","xaxis",xaxis)]]) 
  ny <- max(sapply(xy, length))
  if (length(current_y) == ny) ny <- ny + 1
  y_domains <- lapply(current_y, function(x) current_layout[[x]]$domain)
  main_y_domain <- current_layout$yaxis$domain
  main_y_size <- diff(main_y_domain)
  y_sizes <- sapply(y_domains, diff)
  new_size = main_y_size * size
  new_y_sizes <- c(y_sizes, new_size)
  buffer_size <-  0.04 / sqrt(ny)
  if (side == "top"){
    if (new_size + buffer_size > 1 - max(sapply(y_domains, max))){
      #need to adjust sizes
      norm_factor <- (sum(new_y_sizes) + buffer_size* ny) / (1 - min(sapply(y_domains, min)))
      new_y_sizes <- new_y_sizes / norm_factor
      y_order <- order(sapply(y_domains,min))
      j = min(sapply(y_domains, min))
      for (i in seq_along(y_order)){
        current_layout[[current_y[y_order[i]]]]$domain = c(j,j + new_y_sizes[y_order[i]])
        j = j + new_y_sizes[y_order[i]] + buffer_size
      }
      current_layout[[gsub("y","yaxis",yaxis)]] = c(new_y_layout, 
                                                    list(anchor = xaxis,
                                                         domain = c(j,j + new_y_sizes[length(new_y_sizes)])))
      #Adjust additional y axes as well
      additional_x_axes = names(current_layout)[intersect(grep("xaxis", names(current_layout)), 
                                                          which(names(current_layout) != gsub("x","xaxis",xaxis)))]
      for (k in additional_x_axes){
        additional_y <- Reduce(intersect,list(which(names(current_layout) %in% xy[[k]]), 
                                              which(names(current_layout) != gsub("y","yaxis",yaxis)),
                                              which(names(current_layout) %ni% xy[[gsub("x","xaxis",xaxis)]])))
        if (length(additional_y) == 0) break
        additional_y_domains <- lapply(additional_y, function(x) current_layout[[x]]$domain)
        additional_y_sizes <- sapply(additional_y_domains, diff) 
        additional_new_y_sizes <- additional_y_sizes / norm_factor
        additional_y_bottom <- which(sapply(additional_y_domains, max) < main_y_domain[1])
        additional_y_top <- which(sapply(additional_y_domains, min) > main_y_domain[2])
        additional_y_top_order <- order(sapply(additional_y_domains[additional_y_top],min))
        j = current_layout[["yaxis"]][["domain"]][2] + buffer_size
        for (i in seq_along(additional_y_top_order)){
          current_layout[[additional_y[additional_y_top[additional_y_top_order[i]]]]]$domain = c(j,
                                                                                                 j + additional_new_y_sizes[additional_y_top[additional_y_top_order[i]]])
          j = j + additional_new_y_sizes[additional_y_top[additional_y_top_order[i]]] + buffer_size
        }
        additional_y_bottom_order <- order(sapply(additional_y_domains[additional_y_bottom],max), decreasing = TRUE)
        j = current_layout[["yaxis"]][["domain"]][1] - buffer_size
        for (i in seq_along(additional_y_bottom_order)){
          current_layout[[additional_y[additional_y_bottom[additional_y_bottom_order[i]]]]]$domain = c(j - additional_new_y_sizes[additional_y_bottom[additional_y_bottom_order[i]]],
                                                                                                       j )
          j = j - additional_new_y_sizes[additional_y_bottom[additional_y_bottom_order[i]]] - buffer_size
        }
      }
    } else{
      #Just add new domain
      current_layout[[gsub("y","yaxis",yaxis)]] = c(new_y_layout, 
                                                    list(anchor = xaxis,
                                                         domain = c(max(sapply(y_domains, max)) + buffer_size, 
                                                                    max(sapply(y_domains, max)) + buffer_size + new_size)))
    }
  } else if (side == "bottom"){
    if (new_size + buffer_size > min(sapply(y_domains, min))){
      #need to adjust sizes
      norm_factor <- (sum(new_y_sizes) + buffer_size* ny)  / (max(sapply(y_domains, max)))
      new_y_sizes <- new_y_sizes / norm_factor
      y_order <- order(sapply(y_domains,min))
      j = new_y_sizes[length(new_y_sizes)] + buffer_size
      for (i in seq_along(y_order)){
        current_layout[[current_y[y_order[i]]]]$domain = c(j,j + new_y_sizes[y_order[i]])
        j = j + new_y_sizes[y_order[i]] + buffer_size
      }
      current_layout[[gsub("y","yaxis",yaxis)]] = c(new_y_layout, 
                                                    list(anchor = xaxis,
                                                         domain = c(0,new_y_sizes[length(new_y_sizes)])))
      #Adjust additional y axes as well
      additional_x_axes = names(current_layout)[intersect(grep("xaxis", names(current_layout)), 
                                                          which(names(current_layout) != gsub("x","xaxis",xaxis)))]
      for (k in additional_x_axes){
        additional_y <- Reduce(intersect,list(which(names(current_layout) %in% xy[[k]]), 
                                              which(names(current_layout) != gsub("y","yaxis",yaxis)),
                                              which(names(current_layout) %ni% xy[[gsub("x","xaxis",xaxis)]])))
        if (length(additional_y) == 0) break      
        additional_y_domains <- lapply(additional_y, function(x) current_layout[[x]]$domain)
        additional_y_sizes <- sapply(additional_y_domains, diff) 
        additional_new_y_sizes <- additional_y_sizes / norm_factor
        additional_y_order <- order(sapply(additional_y_domains,min), decreasing = TRUE)
        additional_y_bottom <- which(sapply(additional_y_domains, max) < main_y_domain[1])
        additional_y_top <- which(sapply(additional_y_domains, min) > main_y_domain[2])
        additional_y_top_order <- order(sapply(additional_y_domains[additional_y_top],min))
        j = current_layout[["yaxis"]][["domain"]][2] + buffer_size
        for (i in seq_along(additional_y_top_order)){
          current_layout[[additional_y[additional_y_top[additional_y_top_order[i]]]]]$domain = c(j,
                                                                                                 j + additional_new_y_sizes[additional_y_top[additional_y_top_order[i]]])
          j = j + additional_new_y_sizes[additional_y_top[additional_y_top_order[i]]] + buffer_size
        }
        additional_y_bottom_order <- order(sapply(additional_y_domains[additional_y_bottom],max), decreasing = TRUE)
        j = current_layout[["yaxis"]][["domain"]][1] - buffer_size
        for (i in seq_along(additional_y_bottom_order)){
          current_layout[[additional_y[additional_y_bottom[additional_y_bottom_order[i]]]]]$domain = c(j - additional_new_y_sizes[additional_y_bottom[additional_y_bottom_order[i]]],
                                                                                                       j )
          j = j - additional_new_y_sizes[additional_y_bottom[additional_y_bottom_order[i]]] - buffer_size
        }
      }
    } else{
      #Just add new domain
      current_layout[[gsub("y","yaxis",yaxis)]] = c(new_y_layout, 
                                                    list(anchor = xaxis,
                                                         domain = c(min(sapply(y_domains, min)) - buffer_size - new_size, 
                                                                    min(sapply(y_domains, min)) - buffer_size)))
      
    }
  }
  return(current_layout)
}





