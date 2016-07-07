# Plot Utils, not exported


no_axis = list(title = "",
               zeroline = FALSE,
               showline = FALSE,
               showticklabels = FALSE,
               showgrid = FALSE,
               ticks = "")

dcolorscale <- function(x = 2){
  
  if (x == 1){
    cols = c("#1B9E77")
  } else if (x == 2){
    cols = c("#1B9E77","#7570B3")
  } else if ( x <= 8){
    cols = RColorBrewer::brewer.pal(x, "Dark2")
  } else{
    cols = rainbow(x)
  }
  
  br = rep(seq(0,1,length.out = x  + 1),each = 2)[2:(2*x + 1)]
  
  out <- data.frame(br, rep(cols, each = 2))
  colnames(out) = NULL
  return(out)
}
