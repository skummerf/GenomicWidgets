
setMethod("show", "GenomeTrackWidget",
          function(object){
            print(to_widget(object))
          })


setMethod("show", "GenomicWidget",
          function(object){
            print(to_widget(object))
          })


#' knit_print.GenomicWidget
#' 
#' @param x GenomicWidget object
#' @param options knitr options
#' @keywords internal
#' @export
#' @importFrom knitr knit_print
knit_print.GenomicWidget <- function(x, options){
  knit_print(to_widget(x), options = options)
}


