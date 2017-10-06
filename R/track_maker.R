# # Helper functions
# 
# select_colors <- function(colors, n){
#   if (is.null(colors)){
#     if (n == 1){
#       colors <- "black"
#     } else if ( n <= 8){
#       colors <- RColorBrewer::brewer.pal(n,"Dark2")
#     } else if (n <= 12){
#       colors <- RColorBrewer::brewer.pal(n,"Paired")
#     } else{
#       colors <- rainbow(n)
#     }
#   }
#   return(colors)
# }
# 
# make_ynames <- function(x, start = 1L){
#   stopifnot(x >= 1)
#   end <- start + x - 1
#   out <- paste0("yaxis", seq(start, end))
#   if (start == 1L) out[1] <- "yaxis"
#   out
# }
# 
# 
# get_single_track_trace <- function(i, view, yax, params){
#   
#   signal <- as.vector(assays(make_coverage_matrix(params@data[i], view@range,
#                              binsize = binsize))[[1]])
#   
#   list(
#     x = 
#       seq(
#         relative_position(view,
#                           start(view@range)), 
#         relative_position(view, 
#                           end(view@range)), 
#         length.out = length(tmp_signal)),
#     y = signal,
#     text = params@track_names[i],
#     name = params@track_names[i],
#     hoverinfo = 'x+y+text',
#     line = list(color = x@color[i]),
#     fill = params@fill,
#     mode = params@mode,
#     legendgroup = params@track_names[i],
#     showlegend = params@showlegend,
#     yaxis = gsub("yaxis","y",yax),
#     xaxis = "xaxis")
#   
# }
# 
# get_group_track_trace <- function(group,  yax, view, params){
#   
#   lapply(which(params@group == group),
#          get_single_track_trace,
#          params = params,
#          view = view,
#          yax = yax)
#   
# }
# 
# get_annotation_trace <- function(view, params, yax){
#   
#   # subset transcript info
#   tx_info <- subset_transcripts(view@range, params@annotation)
#   
#   # add stepping
#   if (length(tx_info) >= 1)
#     tx_info <- add_stepping(tx_info)
#   
#   anno_data <- as.data.frame(tx_info, row.names = NULL)
#   if (nrow(anno_data) == 0) return(NULL)
#   anno_data$text <- paste0("Tx ID: ",
#                            anno_data$transcript,
#                            "<br>",
#                            anno_data$feature,
#                            "<br>",
#                            "strand: ",
#                            anno_data$strand)
#   anno_data$start <- relative_position(view, anno_data$start)
#   anno_data$end <- relative_position(view, anno_data$end)
#   anno_data$midpoint <- (anno_data$start + anno_data$end) / 2
#   
#   base_list <- list(yaxis = gsub("yaxis","y",yax),
#                     xaxis = "xaxis",
#                     hoverinfo = "x+text",
#                     opacity = 0,
#                     type = "scatter",
#                     showlegend = FALSE,
#                     mode = "markers")
#   traces <- lapply(unique(anno_data$transcript), function(tname){
#     ix <- which(anno_data$transcript == tname)
#     c(base_list, list(x = anno_data$midpoint[ix],
#                       y = anno_data$stepping[ix],
#                       name = tname))
#     
#   })
#   traces
#   
# }
# 
# get_single_locus_trace <- function(view, ystart, params){
#   
#   if (is.null(params@annotation)){
#     ynames <- make_ynames(1, ystart)
#     traces <- lapply(seq_along(params@data),
#          get_single_track_trace,
#          params = params,
#          view = view,
#          yax = ynames)
#   } else if (params@annotation == "top"){
#     ynames <- make_ynames(2, ystart)
#     traces <- c(get_annotation_trace(view, params, ynames[1]),
#                 lapply(seq_along(params@data),
#                        get_single_track_trace,
#                        params = params,
#                        view = view,
#                        yax = ynames[2]))
#     
#   } else{
#     ynames <- make_ynames(2, ystart)
#     traces <- c(lapply(seq_along(params@data),
#                        get_single_track_trace,
#                        params = params,
#                        view = view,
#                        yax = ynames[1]),
#                 get_annotation_trace(view, params, ynames[2]))
#   }
#   traces
# }
# 
# get_multi_locus_trace <- function(params, views){
#   
#   if (length(views) > 1){
#     
#     if (is.nulll(params@annotation)){
#       ystarts <- seq_along(views)
#     } else{
#       ystarts <- seq(1, length(views) * 2, 2)
#     }
#     
#     out <- mapply(get_single_locus_trace, view = views, ystart = ystarts, 
#          MoreArgs = list(params = params))
#     
#   } else{
#   
#     if (is.null(params@annotation)){
#       ynames <- make_ynames(nlevels(params@groups), ystart)
#       out <- do.call(c,
#                      mapply(get_group_track_trace,
#                             group = levels(params@groups),
#                             yax = ynames, 
#                             MoreArgs = list(params = params,
#                                             view = view))
#       )
#     } else if (params@annotation == "top"){
#       ynames <- make_ynames(nlevels(params@groups)+1, ystart)
#       out <- c(get_annotation_trace(view, params, ynames[1]),
#                do.call(c,
#                        mapply(get_group_track_trace,
#                               group = levels(params@groups),
#                               yax = ynames[seq_len(nlevels(params@groups)) + 1], 
#                               MoreArgs = list(params = params,
#                                               view = view))))
#       
#     } else{
#       ynames <- make_ynames(nlevels(params@groups)+1, ystart)
#       out <- c(do.call(c,
#                        mapply(get_group_track_trace,
#                               group = levels(params@groups),
#                               yax = ynames[seq_len(nlevels(params@groups))], 
#                               MoreArgs = list(params = params,
#                                               view = view))),
#                get_annotation_trace(view, params, ynames[length(ynames)]))
#     }
#   }
#   out
# }
# 
# 
# track_widget <- function(data,
#                          annotation = NULL,
#                          track_names = if (!is.null(names(data)))
#                                            names(data) 
#                                            else basename(data),
#                          groups = track_names,
#                          share_y = TRUE,
#                          showlegend = TRUE,
#                          colors = NULL,
#                          fill = c("tozeroy","none"),
#                          mode = "lines",
#                          annotation_position = c("bottom","top"),
#                          annotation_size = 0.25,
#                          summary = NULL,
#                          layout = list()){
#   
#   params <- list(Class = "TrackParameters",
#                  data = data,
#                  annotation = unpack_transcripts(annotation),
#                  track_names = track_names,
#                  groups = as.factor(groups),
#                  summary = summary,
#                  showlegend = showlegend,
#                  share_y = share_y,
#                  colors = select_colors(colors,length(data)),
#                  fill = match.arg(fill),
#                  annotation_position = match.arg(annotation_position),
#                  mode = mode,
#                  annotation_size = annotation_size,
#                  layout = layout)
#   
#   do.call(new, params)
#   
# }
# 
# summary_y_layout <- function(){
#   
#   
# }
# 
# locus_y_layout <- function(i, ){
#   
#   
# }
# 
# tracks_ylayout <- function(yname, domain, title, range){
#   
#   out <- list()
#   
#   # y axis settings
#   out[[yname]] <- list(title = title,
#                        domain = domain,
#                        zeroline = FALSE,
#                        showgrid = FALSE,
#                        range = range)
#   
#   return(out)
# }
# 
# annotation_ylayout <- function(yname, domain){
#   
#   out <- list()
#   
#   # y axis settings
#   out[[yname]] <- list(domain = domain,
#                        zeroline = FALSE,
#                        showline = FALSE,
#                        showticklabels = FALSE,
#                        showgrid = FALSE,
#                        ticks = "")
#   
#   return(out)
# }
# 
# setup_layout <- function(params, windows, offset, locus_names){
#   
#   out <- list()
#   # Y axes
#   
#   if (length(views) == 1){
#     
#     if (is.null(params@annotations)){
#       ynames <- make_ynames(nlevels(params@groups), 1)
#       domain_end <- 1
#       i <- 1
#       ysize <- 1 / nlevels(params@groups)
#       for (group in levels(params@groups)){
#         out[[ynames[i]]] <- list(title = group,
#                                  domain = c(domain_end - ysize * 0.95,
#                                             domain_end),
#                                  zeroline = FALSE,
#                                  showgrid = FALSE,
#                                  range = c(start(views), end(views)))
#         i <- i + 1
#         domain_end <- domain_end - ysize
#         
#       }
#     } else {
#        if (params@annotation_position == "top"){
#          domain_end <- 1
#          i <- 1
#          ysize <- 1 / length(views)
#          for (view in views){
#            out[[ynames[i]]] <- list(title = locus_names[i],
#                                     domain = c(domain_end - ysize * 0.95,
#                                                domain_end),
#                                     zeroline = FALSE,
#                                     showgrid = FALSE,
#                                     range = c(start(views), end(views)))
#            i <- i + 1
#            domain_end <- domain_end - ysize * params@annotation_size
#            out[[ynames[i]]] <- list(title = locus_names[i],
#                                     domain = c(domain_end - ysize * 0.95,
#                                                domain_end),
#                                     zeroline = FALSE,
#                                     showgrid = FALSE,
#                                     range = c(start(views), end(views)))
#            i <- i + 1
#            domain_end <- domain_end - ysize * (1 - params@annotation_size)
#          }
#          
#          
#        } else{
#          
#        }
#     }
#     
#   } else{
#     if (is.null(params@annotations)){
#       
#     } else {
#       if (params@annotation_position == "top"){
#         
#       } else{
#         
#       }
#     }
#     
#   }
#   
#   
#   if (is.null(params@summary)){
#     # Add xaxis
#     out[["xaxis"]] <- list(zeroline = FALSE,
#                           anchor = gsub("yaxis","y",
#                                         ynames[length(ynames)]),
#                           domain = c(0,1))
#       
#       
#   } else{
#     out[["xaxis"]] <- list(zeroline = FALSE,
#                            anchor = gsub("yaxis","y",tail(ynames,1)),
#                            domain = c(0,0.95*(1 - params@summary@width)))
#     
#     out[["xaxis2"]] <- list(zeroline = FALSE,
#                            anchor = gsub("yaxis","y",tail(ynames2,1)),
#                            domain = c(1 - params@summary@width,1))
#    
#      
#   }
#   
# }
# 
# 
# plot_track <- function(window,
#                        params){
#   
#   
#   # Traces
#   traces <- get_multi_locus_trace(params, views)
#   
#   if (!is.null(summary)){
#     summary_traces <- list() #FILL IN
#     traces <- c(traces, summary_traces)
#   }
#   
#   # Layout
#   
#   layout <- setup_layout(params, views)
#   
#   # Shapes
#   
#   if (!is.null(annotation)){
#     layout$shapes <- make_annotation_shapes(params, views)
#   }
#   
#   
# }


