
          
setMethod(single_locus_view,
          c("GRanges","character"),
          function(window, object, annotation = NULL, 
                   ..., 
                   track_names = ifelse(!is.null(names(object)),
                                        names(object),
                                        basename(object)),
                   groups = NULL,
                   share_y = FALSE,
                   fill = c('tozeroy','none'), 
                   relative = FALSE, 
                   offset = width(window) %/% 2,
                   showlegend = TRUE, 
                   colors = NULL, 
                   mode = 'lines',
                   annotation_position = c("top","bottom"),
                   annotation_size = 0.5){

            
            track_maker <- purrr::partial(make_signal_track, 
                                          window = window,
                                          fill = fill, 
                                          mode = mode, 
                                          showlegend = showlegend
            )
            
            if (!is.null(groups)){
              object_grouped <- split(object, groups)
              track_names_grouped <- split(track_names, groups)
              
              sm <- length(object)
              if (is.null(colors)){
                if (sm == 1){
                  colors <- "black"
                } else if (sm <= 8){
                  colors <- RColorBrewer::brewer.pal(sm,"Dark2")
                } else if (sm <= 12){
                  colors <- RColorBrewer::brewer.pal(sm,"Paired")
                } else{
                  colors <- rainbow(sm)
                }
              }
              
              colors_grouped <- lapply(split(seq_along(colors),groups),
                                       function(x){
                                         colors[x]
                                       })
              
              tracks <- purrr::pmap(list(object = object_grouped,
                                         track_names = track_names_grouped,
                                         colors = colors_grouped,
                                         name = names(object_grouped)), 
                                    track_maker)
            } else{
              if (is.null(colors)){ colors <- rep("black",length(object))}
              
              tracks <- purrr::pmap(list(object = object,
                                         track_names = track_names,
                                         colors = colors), 
                                    track_maker)
              
            }
            
            if (!is.null(annotation)){
              if (match.arg(annotation_position) == "top"){
                heights <- c(annotation_size, rep(1, length(tracks)))
                tracks <- c(make_annotation_track(window,annotation),
                           unname(tracks))
              } else{
                heights <- c(rep(1, length(tracks)), annotation_size)
                tracks <- c(unname(tracks),make_annotation_track(window,
                                                                 annotation))
              }
            } else{
              heights <- rep(1, length(tracks))
            }
            
            if (relative){
             
              if (as.character(strand(window)) == "-"){
                reference <- end(window) - offset
              } else{
                reference <- start(window) + offset
              }
                
              view <- new("ViewRange", range = window, relative = TRUE, 
                         reference = reference)
            } else{
              view <- new("ViewRange", range = window, relative = FALSE)
            }
            
            out <- new("LocusView", as(tracks,"SimpleList"), share_y = share_y,
                         view = view, heights = heights)
            
            return(out)  
            
          })





yaxis_names <- function(x, start = 1L){
 stopifnot(length(x) >= 1)
 end <- start + length(x) - 1
 out <- paste0("yaxis", seq(start, end))
 if (start == 1L) out[1] <- "yaxis"
 out
}

setMethod(get_layout, "AnnotationPlot",
          function(object, yname, domain, ...){
            
            out <- list()
            
            # y axis settings
            out[[yname]] <- list(title = object@trackname,
                                domain = domain,
                                zeroline = FALSE,
                                showline = FALSE,
                                showticklabels = FALSE,
                                showgrid = FALSE,
                                ticks = "")

            return(out)
          })

setMethod(get_layout, "SignalPlot",
          function(object, yname, domain, range, ...){
            
            out <- list()
            
            # y axis settings
            out[[yname]] <- list(title = object@trackname,
                                domain = domain,
                                zeroline = FALSE,
                                showgrid = FALSE,
                                range = range)
            
            return(out)
            
            
          })




setMethod(make_trace, signature = c(x = "LocusView"),
          definition = function(x, ynames, ...){
            traces <- unlist(purrr::map2(as.list(x), ynames, make_trace, 
                                         view = x@view), 
                             recursive = FALSE)
            
            traces
          })

setMethod(make_shapes, signature = c(x = "LocusView"),
          definition = function(x, ynames, ...){
            shapes <- unlist(purrr::map2(as.list(x), ynames, make_shapes, 
                                         view = x@view), 
                             recursive = FALSE)
            
            shapes
          })

setMethod(get_layout, signature = c(object = "LocusView"),
          definition = function(object, ynames, domains, range = NULL, ...){
            unlist(purrr::pmap(list(object = as.list(object), 
                                    yname = ynames,
                                    domain = domains),
                               get_layout,
                               range = range), recursive = FALSE)
          })


