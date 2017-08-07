setMethod(multi_locus_view,
          c("GRange","character"),
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
            new("MultiLocusView", as(single_views,"SimpleList"), share_y = share_y)
            
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


multi_locus_to_plotly_list <- function(x){
  
  lengths <- sapply(x, length)
  ynames <- purrr::map2(as.list(x),
                        cumsum(lengths) - lengths[1] + 1,
                        yaxis_names)
  ynames_flat <- unlist(ynames)
  
  traces <- unlist(purrr::map2(as.list(x),
                        ynames,
                        make_trace), recursive = FALSE)
  
  if (x@share_y){
    range = c(min(x), max(x))
  } else{
    range = NULL
  }
  layout_setting <- list(xaxis = 
                           list(title = "Relative Position",
                                zeroline = FALSE,
                                #showline = FALSE,
                                anchor = gsub("yaxis","y",ynames_flat[length(ynames_flat)]),
                                range = c(relative_position(x[[1]]@view, start(x[[1]]@view)),
                                          relative_position(x[[1]]@view, end(x[[1]]@view)))))
  
  sizes = unlist(purrr::map(as.list(x), function(y) y@heights ))
  
  sizes = sizes / sum(sizes)
  
  domains = list()
  start_domain <- 0
  k <- length(sizes)
  for (i in rev(seq_along(x))){
    domains[[i]] =  list()
    for (j in rev(seq_along(x[[i]]))){
      domains[[i]][[j]] <- c(start_domain, start_domain + (sizes[k]*0.95))
      start_domain <- start_domain + sizes[k]
      k <- k - 1
    }
  }
  
  layout_setting <- c(layout_setting, 
                      unlist(purrr::pmap(list(object = as.list(x), 
                                              yname = ynames,
                                              domain = domains),
                                         get_layout,
                                         range = range), recursive = FALSE))
  
  shapes <- unlist(purrr::map2(as.list(x), ynames, make_shapes), 
                   recursive = FALSE)
  
  layout_setting$shapes <- shapes
  
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

