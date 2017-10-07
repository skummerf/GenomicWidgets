select_colors <- function(colors, n){
  if (is.null(colors)){
    if (n == 1){
      colors <- "black"
    } else if ( n <= 8){
      colors <- RColorBrewer::brewer.pal(n,"Dark2")
    } else if (n <= 12){
      colors <- RColorBrewer::brewer.pal(n,"Paired")
    } else{
      colors <- rainbow(n)
    }
  }
  return(colors)
}


#' max, min for GenomicWidgets objects
#' 
#' @param x object
#' @param na.rm remove na
#' @param ... additional arguments
#' @return numeric
#' @aliases min,AnnotationPlot-method max,AnnotationPlot-method
#' min,SignalPlot-method max,SignalPlot-method min,LocusView-method 
#' max,LocusView-method min,LocusViewList-method max,LocusViewList-method
#' @keywords internal
setMethod(min, "AnnotationPlot",
          function(x, ...){
            NA
          })

setMethod(min, "SignalPlot",
          function(x, na.rm = TRUE){
            
            min(vapply(x@signal, min, 0, na.rm = na.rm))
            
          })


setMethod(max, "SignalPlot",
          function(x, na.rm = TRUE){
            max(vapply(x@signal, max, 0, na.rm = na.rm))
          })

setMethod(max, "AnnotationPlot",
          function(x, ...){
            NA
          })

setMethod(max, "LocusView",
          function(x){
            max(vapply(x, max, 0, na.rm = TRUE), na.rm = TRUE)
          })

setMethod(min, "LocusView",
          function(x){
            min(vapply(x, min, 0, na.rm = TRUE), na.rm = TRUE)
          })

setMethod(max, "LocusViewList",
          function(x){
            max(vapply(x, max, 0, na.rm = TRUE), na.rm = TRUE)
          })

setMethod(min, "LocusViewList",
          function(x){
            min(vapply(x, min, 0, na.rm = TRUE), na.rm = TRUE)
          })

