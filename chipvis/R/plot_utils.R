# Some utility functions for use with plotly
#' @importFrom RColorBrewer brewer.pal brewer.pal.info

no_axis = list(title = "",
               zeroline = FALSE,
               showline = FALSE,
               showticklabels = FALSE,
               showgrid = FALSE,
               ticks = "")


brewer.pal.helper <- function(n, name){
  if (n < 3){
    return(RColorBrewer::brewer.pal(3,name)[seq_len(n)])
  } else if (n > RColorBrewer::brewer.pal.info[name, "maxcolors"]){
    return(rep(RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[name, "maxcolors"],name),
               length.out = n))
  } else {
    return(RColorBrewer::brewer.pal(n,name))
  }
}


# makes a discrete color scale for plotly
# palette should be an RColorBrewer palette name
# x is number of colors desired

#' dcolorscale
#' 
#' @param x number of items
#' @param palette name of RColorBrewer palette
#' @param cols vector of color names, overrides palette if given
#' @return data.frame of colors and breaks to give as colorscale input to plotly
#' @export
#' @author Alicia Schep
dcolorscale <- function(x = 2, palette = "Dark2", cols = NULL){
  
  stopifnot(x >= 1)
  
  if (is.null(cols)){
    if (x == 1){
      cols = RColorBrewer::brewer.pal(3, palette)[1]
    } else if (x == 2){
      cols = RColorBrewer::brewer.pal(3, palette)[c(1,3)]
    } else if ( x <= RColorBrewer::brewer.pal.info[palette, "maxcolors"]){
      cols = RColorBrewer::brewer.pal(x, palette)
    } else{
      warning("Desired pallete contains insufficient colors")
      cols = rainbow(x)
    }    
  } else{
    cols = cols[1:x]
  }

  br = rep(seq(0,1,length.out = x  + 1),each = 2)[2:(2*x + 1)]
  
  out <- data.frame(br, rep(cols, each = 2))
  colnames(out) = NULL
  return(out)
}

alternating_colorscale <- function(x = 2, cols = c("lightgray","darkgray")){
  
  stopifnot(x >= 1)
  
  cols = rep(cols, (x + 1) %/% 2)
  cols = cols[1:x]
  
  br = rep(seq(0,1,length.out = x  + 1),each = 2)[2:(2*x + 1)]
  
  out <- data.frame(br, rep(cols, each = 2))
  colnames(out) = NULL
  return(out)
}

choose_discrete_palette <- function(x, existing = c()){
  
  qual_colors <- RColorBrewer::brewer.pal.info[which(RColorBrewer::brewer.pal.info$category == "qual"),]
  new <- which(rownames(qual_colors) %ni% existing)
  if (length(new) == 0){
    existing_factors <- as.factor(existing)
    existing_tab <- tabulate(existing_factors)
    new <- which(rownames(qual_colors) %in% levels(existing_factors)[which(existing_tab < max(existing_tab))])
  }
  enough <- which(qual_colors$maxcolors >= x)
  if (length(intersect(new,enough)) > 0){
    qual_colors <- qual_colors[intersect(new,enough),]
    cbs <- which(qual_colors$colorblind)
    if (length(cbs) > 0){
      qual_colors <- qual_colors[cbs,]
    }
  } else if (length(enough) > 0){
    qual_colors <- qual_colors[enough,]
    cbs <- which(qual_colors$colorblind)
    if (length(cbs) > 0){
      qual_colors <- qual_colors[cbs,]
    }
  } else if (length(new) > 0){
    qual_colors <- qual_colors[new,]
    cbs <- which(qual_colors$colorblind)
    if (length(cbs) > 0){
      qual_colors <- qual_colors[cbs,]
    }
  } else{
    cbs <- which(qual_colors$colorblind)
    if (length(cbs) > 0){
      qual_colors <- qual_colors[cbs,]
    }
  }
  return(rownames(qual_colors)[1])  
}


# ccolorscale <- function(z, palette = "Reds", mid = NULL, zmin = min(z), zmax = max(z)){
#   vals <- unique(scales::rescale(c(z)))
#   if (!is.null(mid)){
#     if (sym){
#       domain <- c(min(zmin, mid - (zmax - mid)),
#                   max(zmaxz, mid + zmaxz))
#       vals2 <- unique(scales::rescale_mid(z, from = domain,
#                                           mid = mid))
#     } else {
#       vals2 <- unique(scales::rescale_mid(z, mid = mid))
#     }
#   } else {
#     domain <- min()
#     vals2 <- vals
#   }
#   o <- order(vals, decreasing = FALSE)
#   cols <- scales::col_numeric(palette)(vals2)
#   colz <- setNames(data.frame(vals[o], cols[o]), NULL) 
#  return(colz)
# }


# choose_continuous_palette(signal, existing){
#   if (min(scale) < 0 && max(scale) > 0){
#     cont_colors <- RColorBrewer::brewer.pal.info[which(RColorBrewer::brewer.pal.info$category == "seq"),]
#   }
#   
#   
#   new <- which(rownames(cont_colors) %ni% existing)
#   if (length(new) == 0){
#     existing_factors <- as.factor(existing)
#     existing_tab <- tabulate(existing_factors)
#     new <- which(rownames(cont_colors) %in% levels(existing_factors)[which(existing_tab < max(existing_tab))])
#   }
#   enough <- which(cont_colors$maxcolors >= x)
#   if (length(intersect(new,enough)) > 0){
#     cont_colors <- cont_colors[intersect(new,enough),]
#     cbs <- which(cont_colors$colorblind)
#     if (length(cbs) > 0){
#       cont_colors <- cont_colors[cbs,]
#     }
#   } else if (length(enough) > 0){
#     cont_colors <- cont_colors[enough,]
#     cbs <- which(cont_colors$colorblind)
#     if (length(cbs) > 0){
#       cont_colors <- cont_colors[cbs,]
#     }
#   } else if (length(new) > 0){
#     cont_colors <- cont_colors[new,]
#     cbs <- which(cont_colors$colorblind)
#     if (length(cbs) > 0){
#       cont_colors <- cont_colors[cbs,]
#     }
#   } else{
#     cbs <- which(cont_colors$colorblind)
#     if (length(cbs) > 0){
#       cont_colors <- cont_colors[cbs,]
#     }
#   }
#   return(rownames(cont_colors)[1])  
# }


# makes x based on colnames of mat if available
# if not available, just uses 1 to number of columns
default_x <- function(mat){
  if (is.null(colnames(mat))){
    return(1:ncol(mat))
  } else{
    colnames(mat)
  }
}
# makes y based on rownames of mat if available
# if not available, just uses 1 to number of rows
default_y <- function(mat){
  if (is.null(rownames(mat))){
    return(1:nrow(mat))
  } else{
    rownames(mat)
  }
}
