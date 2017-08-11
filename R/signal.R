

setMethod("make_signal_track", c("GRanges","character"),
          function(window, object, bin.num = 1000, ..., 
                   track_names = ifelse(!is.null(names(object)),
                                        basename(names(object)),
                                        object),
                   fill = c('tozeroy','none'), 
                   showlegend = TRUE, 
                   colors = NULL, 
                   mode = 'lines',
                   name = ifelse(length(track_names) > 1, "Coverage", track_names)){
            
            make_signal_track(as(window,"ViewRange"), object, bin_num = bin_num,
                              ..., 
                              track_names = track_names, 
                              fill = match.arg(fill), showlegend = showlegend,
                              colors = colors, mode = mode, name = name)
          })

setMethod("make_signal_track", c("ViewRange","character"),
          function(window, object, bin.num = 1000, ..., 
                   track_names = ifelse(!is.null(names(object)),
                                        basename(names(object)),
                                        object),
                   fill = c('tozeroy','none'), 
                   showlegend = TRUE, 
                   colors = NULL, 
                   mode = 'lines',
                   name = "Coverage"){
            
            fill <- match.arg(fill)
            
            sm <- genomation::ScoreMatrixList(object, as(window,"GRanges"), 
                                              bin.num = bin.num, ...)
            names(sm) <- track_names
            
            
            if (is.null(colors)){
              if (length(sm == 1)){
                colors = "black"
              } else if (length(sm <= 8)){
                colors = RColorBrewer::brewer.pal(length(sm),"Dark2")
              } else if (length(sm <= 12)){
                colors = RColorBrewer::brewer.pal(length(sm),"Paired")
              } else{
                colors = rainbow(length(sm))
              }
            }
            
            # Make SignalPlot object
            new("SignalPlot",
                signal = sm,
                color = colors,
                mode = mode,
                fill = fill,
                trackname = name,
                showlegend = showlegend)
          })

#' @export
setMethod(make_trace, signature = c(x = "SignalPlot"),
          definition = function(x, yax, view, xax = "xaxis", ...){
            
            trace_data <- lapply(seq_len(length(x@signal)),
                                 function(i){
                                   tmp_signal = as.vector(x@signal[[i]])
                                   list(x = seq(relative_position(view,
                                                                  start(view)), 
                                                relative_position(view, 
                                                                  end(view)), 
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

signal_to_plotly_list <- function(x){
  traces <- make_trace(x, "yaxis")
  layout_setting <- list()#get_layout(x)
  out <- list(data = traces,
              layout = layout_setting,
              source = "Annotation Track",#,x@source,
              config = list(modeBarButtonsToRemove =
                              c("sendDataToCloud",
                                "autoScale2d")))
  attr(out, "TOJSON_FUNC") <- function(x, ...) {
    jsonlite::toJSON(x, digits = 50, auto_unbox = TRUE, force = TRUE,
                     null = "null", na = "null", ...)
  }
  out
}

#' @export
setMethod(to_widget,
          c("SignalPlot"),
          function(p){
            out <- signal_to_plotly_list(p)
            htmlwidgets::createWidget(
              name = "chipVis",
              x = out,
              width = out$layout$width,
              height = out$layout$height,
              sizingPolicy = htmlwidgets::sizingPolicy(browser.fill = TRUE,
                                                       viewer.fill = TRUE,
                                                       defaultWidth = "100%",
                                                       defaultHeight = 400),
              dependencies = plotly_dependency())
          })





