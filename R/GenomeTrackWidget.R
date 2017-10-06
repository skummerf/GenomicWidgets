

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



#' knit_print.GenomeTrackWidget
#' 
#' @param x GenomeTrackWidget object
#' @param options knitr options
#' @keywords internal
#' @export
#' @importFrom knitr knit_print
knit_print.GenomeTrackWidget <- function(x, options){
  knit_print(to_widget(x), options = options)
}



#' knit_print.LocusViewList
#' 
#' @param x LocusViewList object
#' @param options knitr options
#' @keywords internal
#' @export
#' @importFrom knitr knit_print
knit_print.LocusViewList <- function(x, options){
  knit_print(to_widget(x), options = options)
}

#' knit_print.LocusSummaryList
#' 
#' @param x LocusSummaryList object
#' @param options knitr options
#' @keywords internal
#' @export
#' @importFrom knitr knit_print
knit_print.LocusSummaryList <- function(x, options){
  knit_print(to_widget(x), options = options)
}

