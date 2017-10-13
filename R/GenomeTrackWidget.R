

plotly_dependency <- function(){
  htmltools::htmlDependency(
    "plotlyjs", "1.29.2",
    src = system.file('htmlwidgets', 'lib', 'plotlyjs', package = 'iheatmapr'),
    script = "plotly-latest.min.js",
    stylesheet = "plotly-htmlwidgets.css"
  )
}

setMethod("show", "GenomeTrackWidget",
          function(object){
            print(to_widget(object))
          })


setMethod("show", "LocusViewList",
          function(object){
            print(to_widget(object))
          })

setMethod("show", "LocusSummaryList",
          function(object){
            print(to_widget(object))
          })



#' knit_print for GenomicWidgets objects
#' 
#' knit_print method for objects that generate GenomicWidgets
#' @param x GenomeTrackWidget object
#' @param options knitr options
#' @return htmlwidgets object
#' @keywords internal
#' @export
#' @rdname knit_print_GenomicWidgets
#' @name knit_print_GenomeWidgets
#' @importFrom knitr knit_print
knit_print.GenomeTrackWidget <- function(x, options){
  knit_print(to_widget(x), options = options)
}



#' @rdname knit_print_GenomicWidgets
#' @name knit_print_GenomeWidgets#' 
#' @export
knit_print.LocusViewList <- function(x, options){
  knit_print(to_widget(x), options = options)
}

#' @rdname knit_print_GenomicWidgets
#' @name knit_print_GenomeWidgets#' 
#' @export
knit_print.LocusSummaryList <- function(x, options){
  knit_print(to_widget(x), options = options)
}

