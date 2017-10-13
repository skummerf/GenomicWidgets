#' Shiny bindings for GenomicWidgets
#' 
#' Output and render functions for using GenomicWidgets within Shiny 
#' 
#' @param outputId output variable to read from
#' @param width,height Must be a valid CSS unit (like \code{"100\%"},
#'   \code{"400px"}, \code{"auto"}) or a number, which will be coerced to a
#'   string and have \code{"px"} appended.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This 
#'   is useful if you want to save an expression in a variable.
#' @return analagous to \code{\link[shiny]{renderPlot}} and 
#'   \code{\link[shiny]{plotOutput}} from shiny package, renderGenomicWidgets 
#'   is used in server side code to render the plot, and GenomicWidgetsOutput
#'   is used within the ui code.
#' @importFrom htmlwidgets shinyWidgetOutput shinyRenderWidget
#' @name GenomicWidgets-shiny
#'
#' @export
GenomicWidgetsOutput <- function(outputId, width = "100%", height = "400px") {
  htmlwidgets::shinyWidgetOutput(outputId, "GenomicWidgets", width, height, 
                                 package = "GenomicWidgets")
}

#' @param expr An expression that generates an Iheatmap object
#' @param env The environment in which to evaluate \code{expr}.
#' @rdname GenomicWidgets-shiny
#' @export
renderGenomicWidgets <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  expr <- call("to_widget", expr)
  htmlwidgets::shinyRenderWidget(expr, GenomicWidgetsOutput, env, quoted = TRUE)
}
