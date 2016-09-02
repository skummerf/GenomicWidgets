
#' aggregate_profile_plot
#' 
#' @param mats coverage matrices, as created by \code{\link{make_coverage_matrix}}
#' @param groups vector of groups to which matrices belong
#' @param colors either RColorbrewer pallete name or vector of colors
#' @param positions positions corresponding to each column of one of the matrices
#' @param ylab axis name for y axis
#' @param xlab axis name for x axis
#' @param scale_method method to use for scaling, see  \code{\link{normalize_coverage_matrix}}
#' @param pct percentile argument if PercentileMax chosen as scale_method
#' @param scale_factors scale_factors if "scalar" chosen as scale_method
#' @param showlegend show the legend? 
#' @param source source name, for use in shiny apps
#'
#' @author Alicia Schep
#' @export
aggregate_profile_plot <- function(mats, 
                                   groups = NULL,
                                   colors = "Dark2",
                                   positions = as.numeric(colnames(mats[[1]])),
                                   ylab = "ChIP Signal", 
                                   xlab = "Position",
                                   scale_method = c("none",
                                                    "localRms", 
                                                    "localMean", 
                                                    "localNonZeroMean", 
                                                    "PercentileMax", 
                                                    "scalar"),
                                   pct = 0.95,
                                   scale_factors = 1,
                                   showlegend = TRUE,
                                   source = AGGREGATE_SOURCE){
  
  scale_method = match.arg(scale_method)
  if (scale_method != "none"){
    mats <- normalize_coverage_matrix(mats,
                                      method = scale_method, 
                                      pct = pct, 
                                      scalar = scale_factors)
  } 

  col_aggs <- lapply(mats, colMeans, na.rm = TRUE)
  
  if (!is.null(groups)){
    col = rep(groups, each = length(positions))
  } else{
    col = rep(names(mats), each = length(positions))
  }
  
  df <- data.frame(x = rep(positions, length(col_aggs)),
                   y =   unlist(col_aggs, use.names = FALSE),
                   names = rep(names(mats), each = length(positions)),
                   col = col)
  
  df %>% group_by(names) %>% plot_ly(x = ~x,
          y = ~y,
          color = ~col,
          colors = colors,
          text = ~names,
          type = "scatter",
          mode = "lines",
          source = source) %>% layout(hovermode = "closest",
                                     xaxis = list(title = xlab),
                                     yaxis = list(title = ylab),
                                     showlegend = showlegend)
  
  
}