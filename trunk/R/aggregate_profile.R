
#' aggregate_profile_plot
#' 
#' @author Alicia Schep
#' @export
aggregate_profile_plot <- function(cvg_mats, 
                                   positions = colnames(cvg_mats[[1]]),
                                   ylab = "ChIP Signal", 
                                   xlab = "Position",
                                   showlegend = TRUE){
  
  col_aggs <- lapply(cvg_mats, colMeans)
  
  plot_ly(x = rep(positions, length(col_aggs)),
          y = unlist(col_aggs, use.names = FALSE),
          color = rep(names(col_aggs), each = length(positions)),
          text = rep(names(col_aggs), each = length(positions)),
          mode = "lines",
          source = "agg") %>% layout(hovermode = "closest",
                                     xaxis = list(title = xlab),
                                     yaxis = list(title = ylab),
                                     showlegend = showlegend)
  
  
}