# Some utility functions for use with plotly

no_axis = list(title = "",
               zeroline = FALSE,
               showline = FALSE,
               showticklabels = FALSE,
               showgrid = FALSE,
               ticks = "")

# makes a discrete color scale for plotly
# palette should be an RColorBrewer palette name
# x is number of colors desired

#' dcolorscale
#' 
#' @param x number of items
#' @param palette name of RColorBrewer palette
#' @export
#' @author Alicia Schep
dcolorscale <- function(x = 2, palette = "Dark2"){
  
  if (x == 1){
    cols = RColorBrewer::brewer.pal(3, palette)[1]
  } else if (x == 2){
    cols = RColorBrewer::brewer.pal(3, palette)[c(1,3)]
  } else if ( x <= 8){
    cols = RColorBrewer::brewer.pal(x, palette)
  } else{
    cols = rainbow(x)
  }
  
  br = rep(seq(0,1,length.out = x  + 1),each = 2)[2:(2*x + 1)]
  
  out <- data.frame(br, rep(cols, each = 2))
  colnames(out) = NULL
  return(out)
}

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
