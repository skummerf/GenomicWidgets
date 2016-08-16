
#' aggregate_profile_plot
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
                                   scale_factor = 1,
                                   showlegend = TRUE){
  
  scale_method = match.arg(scale_method)
  if (scale_method != "none"){
    mats <- normalize_coverage_matrix(mats,
                                      method = scale_method, 
                                      pct = pct, 
                                      scalar = scale_factor)
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
          source = "agg") %>% layout(hovermode = "closest",
                                     xaxis = list(title = xlab),
                                     yaxis = list(title = ylab),
                                     showlegend = showlegend)
  
  
}