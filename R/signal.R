# Helper methods for making the signal tracks
# Not exported

setMethod("make_signal_track", c("GRanges","character"),
          function(window, object, binsize = 25, ..., 
                   track_names = ifelse(!is.null(names(object)),
                                        basename(names(object)),
                                        object),
                   fill = c('tozeroy','none'), 
                   showlegend = TRUE, 
                   colors = NULL, 
                   mode = 'lines',
                   name = if (length(track_names) > 1) "Coverage" else 
                     track_names){
            
            fill <- match.arg(fill)
            
            names(object) <- track_names
            sm <- make_coverage_matrix(object, as(window,"GRanges"),
                                       binsize = binsize, ...)

            
            if (is.null(colors)){
              if (length(sm == 1)){
                colors <- "black"
              } else if (length(sm <= 8)){
                colors <- RColorBrewer::brewer.pal(length(sm),"Dark2")
              } else if (length(sm <= 12)){
                colors <- RColorBrewer::brewer.pal(length(sm),"Paired")
              } else{
                colors <- rainbow(length(sm))
              }
            }
            
            # Make SignalPlot object
            new("SignalPlot",
                signal = assays(sm),
                color = colors,
                mode = mode,
                fill = fill,
                trackname = name,
                showlegend = showlegend)
          })


setMethod(make_trace, signature = c(x = "SignalPlot"),
          definition = function(x, yax, view, xax = "xaxis", ...){
            
            trace_data <- lapply(seq_len(length(x@signal)),
                                 function(i){
                                   tmp_signal <- as.vector(x@signal[[i]])
                                   list(
                                     x = 
                                       seq(
                                         relative_position(view,
                                                           start(view@range)), 
                                         relative_position(view, 
                                                           end(view@range)), 
                                         length.out = length(tmp_signal)),
                                     y = tmp_signal,
                                     text = names(x@signal)[i],
                                     name = names(x@signal)[i],
                                     hoverinfo = 'x+y+text',
                                     line = list(color = x@color[i]),
                                     fill = x@fill,
                                     mode = x@mode,
                                     legendgroup = names(x@signal)[i],
                                     showlegend = x@showlegend,
                                     yaxis = gsub("yaxis","y",yax),
                                     xaxis = gsub("xaxis","x",xax))
                                 })
            trace_data
            })

setMethod("make_shapes", c(x = "SignalPlot"),
          function(x, ...){
            return(NULL)
          })







