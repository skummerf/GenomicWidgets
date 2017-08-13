setMethod(multi_locus_view,
          c("GRanges","character"),
          function(windows, 
                   object, 
                   annotation = NULL, 
                   ..., 
                   track_names = ifelse(!is.null(names(object)),
                                        names(object),
                                        basename(object)),
                   name = mcols(windows)$name,
                   share_y = FALSE,
                   fill = c('tozeroy','none'), 
                   showlegend = (length(object) > 1), 
                   colors = NULL, 
                   mode = 'lines',
                   annotation_position = c("bottom","top"),
                   annotation_size = 0.2){
            
            multi_locus_view(as(windows, "RelativeViewRange"),
                               object,
                               annotation,
                               ...,
                               track_names = track_names,
                               name = name,
                               share_y = share_y,
                               fill = match.arg(fill),
                               showlegend = showlegend,
                               colors = colors,
                               mode = mode,
                               annotation_position = match.arg(annotation_position),
                               annotation_size = annotation_size)
          })


setMethod(multi_locus_view,
          c("RelativeViewRange","character"),
          function(windows, 
                   object, 
                   annotation = NULL, ..., 
                   track_names = ifelse(!is.null(names(object)),
                                        names(object),
                                        basename(object)),
                   name = mcols(windows)$name,
                   share_y = FALSE,
                   fill = c('tozeroy','none'), 
                   showlegend = (length(object) > 1), 
                   colors = NULL, 
                   mode = 'lines',
                   annotation_position = c("bottom","top"),
                   annotation_size = 0.2){
            
            annotation_position <- match.arg(annotation_position)
            
            
            sm <- length(object)
            if (is.null(colors)){
              if (sm == 1){
                colors = "black"
              } else if (sm <= 8){
                colors = RColorBrewer::brewer.pal(sm,"Dark2")
              } else if (sm <= 12){
                colors = RColorBrewer::brewer.pal(sm,"Paired")
              } else{
                colors = rainbow(sm)
              }
            }
            
            single_views <- purrr::map(seq_along(windows),
                                       function(x){
                                         single_locus_view(windows[x],
                                                           object = object,
                                                           annotation = annotation,
                                                           track_names = track_names,
                                                           groups = rep(name[x],length(object)),
                                                           fill = fill,
                                                           showlegend = if (x == 1) showlegend else FALSE,
                                                           colors = colors,
                                                           mode = mode,
                                                           annotation_position = annotation_position,
                                                           annotation_size = annotation_size)
                                       })
            ll <- new("LocusViewList", as(single_views,"SimpleList"), share_y = share_y)
            ll
            #new("MultiLocusView", tracks = ll)
          })
          

setMethod(multi_locus_view,
          c("character","character"),
          function(windows, object, annotation = NULL, ..., 
                   track_names = ifelse(!is.null(names(object)),
                                        names(object),
                                        basename(object)),
                   groups = NULL,
                   share_y = FALSE,
                   fill = c('tozeroy','none'), 
                   relative = FALSE, 
                   showlegend = TRUE, 
                   colors = NULL, 
                   mode = 'lines',
                   annotation_position = c("top","bottom"),
                   annotation_size = 0.5){
            
            
            
            
          })      

#' @export
setMethod(make_track_plotter,
          c("character"),
          function(object, 
                   annotation = NULL, 
                   ..., 
                   track_names = ifelse(!is.null(names(object)),
                                        names(object),
                                        basename(object)),
                   groups = NULL,
                   share_y = FALSE,
                   fill = c('tozeroy','none'), 
                   relative = FALSE, 
                   showlegend = TRUE, 
                   colors = NULL, 
                   mode = 'lines',
                   annotation_position = c("bottom","top"),
                   annotation_size = 0.25){
            
            fill = match.arg(fill)
            annotation_position = match.arg(annotation_position)
            
            purrr::partial(multi_locus_view,
                           object = object,
                           annotation = annotation,
                           ...,
                           track_names = track_names,
                           groups = groups,
                           share_y = share_y,
                           fill = fill,
                           relative = relative,
                           showlegend = showlegend,
                           colors = colors,
                           mode = mode,
                           annotation_position = annotation_position,
                           annotation_size = annotation_size)
            
            
          })

setMethod(make_trace, signature = c(x = "LocusViewList"),
          definition = function(x, ynames, ...){
            unlist(purrr::map2(as.list(x),
                               ynames,
                               make_trace), recursive = FALSE)
          })

setMethod(make_trace, signature = c(x = "LocusSummaryList"),
          definition = function(x, ynames, xaxis = "xaxis2", ...){
            unlist(purrr::map2(as.list(x),
                               ynames,
                               make_trace, 
                               xax = xaxis), recursive = FALSE)
          })

setMethod(make_shapes, signature = c(x = "LocusViewList"),
          definition = function(x, ynames, ...){
            unlist(purrr::map2(as.list(x), ynames, make_shapes), 
                   recursive = FALSE)
          })

setMethod(get_layout, signature = c(object = "LocusViewList"),
          definition = function(object, ynames, x_domain = c(0,1), ...){

            ynames_flat <- unlist(ynames)
            
            if (object@share_y){
              range = c(min(object), max(object))
            } else{
              range = NULL
            }
            layout_setting <- list(xaxis = 
                                     list(title = "Relative Position",
                                          zeroline = FALSE,
                                          #showline = FALSE,
                                          anchor = gsub("yaxis","y",ynames_flat[length(ynames_flat)]),
                                          range = c(relative_position(object[[1]]@view, start(object[[1]]@view)),
                                                    relative_position(object[[1]]@view, end(object[[1]]@view))),
                                          domain = x_domain))
            
            sizes = unlist(purrr::map(as.list(object), function(y) y@heights ))
            
            sizes = sizes / sum(sizes)
            
            domains = list()
            start_domain <- 0
            k <- length(sizes)
            for (i in rev(seq_along(object))){
              domains[[i]] =  list()
              for (j in rev(seq_along(object[[i]]))){
                domains[[i]][[j]] <- c(start_domain, start_domain + (sizes[k]*0.95))
                start_domain <- start_domain + sizes[k]
                k <- k - 1
              }
            }
            
            layout_setting <- c(layout_setting, 
                                unlist(purrr::pmap(list(object = as.list(object), 
                                                        yname = ynames,
                                                        domain = domains),
                                                   get_layout,
                                                   range = range), recursive = FALSE))
            
            layout_setting
          })

setMethod(get_layout, signature = c(object = "LocusSummaryList"),
          definition = function(object, ynames, xax = "xaxis2", x_domain = c(0,1), ...){
            

            layout_setting <- list()
            layout_setting[[xax]] = list(zeroline = FALSE,
                                          #showline = FALSE,
                                          anchor = gsub("yaxis","y",ynames[length(ynames)]),
                                         domain = x_domain)
            
            sizes = rep(1 / length(object), length(object))
            domains = list()
            start_domain <- 0
            for (i in rev(seq_along(object))){
              domains[[i]] <- c(start_domain, start_domain + (sizes[i]*0.95))
              start_domain <- start_domain + sizes[i]
            }
            
            layout_setting <- c(layout_setting, 
                                unlist(purrr::pmap(list(object = as.list(object), 
                                                        yname = ynames,
                                                        domain = domains),
                                                   get_layout,
                                                   anchor = xax), 
                                                   recursive = FALSE))
            
            layout_setting
          })




multi_locus_to_plotly_list <- function(x){
  
  if (length(x@tracks) >= 1){
    
    lengths <- sapply(x@tracks, length)
    track_ynames <- purrr::map2(as.list(x@tracks),
                                cumsum(lengths) - lengths[1] + 1,
                                yaxis_names)
    
    traces <- make_trace(x@tracks, track_ynames)
    
    if (length(x@summaries) == 0){
      x_domain = c(0,1)
    } else{
      x_domain = c(0, (1 - x@summary_width) * 0.95)
    }
    
    layout_setting <- get_layout(x@tracks, 
                                 track_ynames,
                                 x_domain = x_domain)
    
    shapes <- make_shapes(x@tracks, track_ynames)
    
    layout_setting$shapes <- shapes
    
    xax <- "xaxis2"
    
  } else{
    lengths <- c()
    traces <- list()
    layout_setting <- list()
    xax <- "xaxis"
  }
  
  if (length(x@summaries) !=0){
    summary_ynames <- yaxis_names(x@summaries, sum(lengths) + 1)
    traces <- c(traces, make_trace(x@summaries, summary_ynames, xax = xax))
    layout_setting <- c(layout_setting, get_layout(x@summaries,
                                                   summary_ynames,
                                                   xax = xax,
                                                   x_domain = c(1 - x@summary_width,
                                                              1)))
    
  }
  
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
          c("LocusViewList"),
          function(p){
            p <- new("MultiLocusView", tracks = p)
            out <- multi_locus_to_plotly_list(p)
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

#' @export
setMethod(to_widget,
          c("LocusSummaryList"),
          function(p){
            p <- new("MultiLocusView", summaries = p, summary_width = 1)
            out <- multi_locus_to_plotly_list(p)
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



#' @export
setMethod(to_widget,
          c("MultiLocusView"),
          function(p){
            out <- multi_locus_to_plotly_list(p)
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

