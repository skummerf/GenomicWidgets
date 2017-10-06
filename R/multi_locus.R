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
                   groups = NULL,
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
                colors <- "black"
              } else if (sm <= 8){
                colors <- RColorBrewer::brewer.pal(sm,"Dark2")
              } else if (sm <= 12){
                colors <- RColorBrewer::brewer.pal(sm,"Paired")
              } else{
                colors <- rainbow(sm)
              }
            }
            
            if (length(windows) == 1){
              
              single_views <- list(single_locus_view(windows,
                                                     object = object,
                                                     annotation = annotation,
                                                     track_names = track_names,
                                                     groups = groups,
                                                     fill = fill,
                                                     showlegend = showlegend ,
                                                     colors = colors,
                                                     mode = mode,
                                                     annotation_position = 
                                                       annotation_position,
                                                     annotation_size = 
                                                       annotation_size))

            } else{

              if (is.null(name)){
                if (!is.null(names(windows))){
                  name <- names(windows)
                } else{
                  name <- seq_along(windows)
                }
              }
              
              single_views <- 
                purrr::map(seq_along(windows),
                           function(x){
                             single_locus_view(
                               windows[x],
                               object = object,
                               annotation = annotation,
                               track_names = track_names,
                               groups = rep(name[x],
                                            length(object)),
                               fill = fill,
                               relative = TRUE,
                               showlegend = 
                                 if (x == 1) showlegend else FALSE,
                               colors = colors,
                               mode = mode,
                               annotation_position = annotation_position,
                               annotation_size = annotation_size)
                           })
              
              
            }
            ll <- new("LocusViewList", as(single_views,"SimpleList"), 
                      share_y = share_y)
            return(ll)
          })





#' make_track_plotter
#' 
#' Function to generate a function that takes in a range and plots coverage track
#' @param object vector of bam or bigwig file names
#' @param annotation TxDb or OrganismDb object
#' @param ... additional arguments
#' @param track_names names to associate with each file
#' @param groups vector of group assignments.  traces will be grouped onto subplots
#' based on group assignments (if only showing 1 region)
#' @param share_y share the y axis?
#' @param fill fillmode for line plot
#' @param relative ignore for now
#' @param showlegend show the legend?
#' @param colors colors for each bam file
#' @param mode mode for plot
#' @param annotation_position plot annotations on bottom or on top of signal traces
#' @param annotation_size relative size of annotation plot
#' @export
#' @rdname make_track_plotter
#' @name make_track_plotter
#' @aliases make_track_plotter,character-method
#' 
#' @author Alicia Schep and Justin Finkle
#' @return function to make interactive track plots
#' @examples 
#' 
#' library(GenomicRanges)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' 
#' ## First we'll read in some sample data
#' genomation_dir <- system.file("extdata", package = "genomationData")
#' samp.file <- file.path(genomation_dir,'SamplesInfo.txt')
#' samp.info <- read.table(samp.file, header=TRUE, sep="\t", 
#'                         stringsAsFactors = FALSE)
#' samp.info$fileName <- file.path(genomation_dir, samp.info$fileName)
#' ctcf.peaks = genomation::readBroadPeak(system.file("extdata",
#'                          "wgEncodeBroadHistoneH1hescCtcfStdPk.broadPeak.gz",
#'                           package = "genomationData"))
#' ctcf.peaks = ctcf.peaks[seqnames(ctcf.peaks) == "chr21"]
#' 
#' ## resize peaks to size 1000
#' ctcf.peaks = resize(ctcf.peaks, width = 10000, fix = "center")
#' 
#' ## Make track plotter
#' 
#' track_plotter <- make_track_plotter(samp.info$fileName[1:3], 
#'   annotation = TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'   track_names = samp.info$sampleName[1:3] , 
#'   share_y = TRUE)
#' if (interactive()){
#'   track_plotter(ctcf.peaks[1])
#'   track_plotter(ctcf.peaks[1:3])
#' }   
#' 
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
            
            fill <- match.arg(fill)
            annotation_position <- match.arg(annotation_position)
            
            default_arglist <- list(
              object = object,
              annotation = unpack_transcripts(annotation),
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
              annotation_size = annotation_size
            )
            
            out <- function(windows, ...){
              arglist <- modifyList(default_arglist, 
                                    list(...))
              do.call(multi_locus_view, c(list(windows = windows), arglist))
              
            }
            
            out
          })

#' Combine genome track view with locus summary
#' 
#' @param track_function function to make tracks, as created by
#'  \code{\link{make_track_plotter}}
#' @param summary_function function to make locus summaries, as created by
#'  \code{\link{make_summary_plotter}}
#' @param windows GenomicRanges of windows that can serve as inputs to 
#' track_function
#' @param row_names vector of row names that can serve as inputs to 
#' summary_function
#' @param summary_width width of summary plot in resulting visualization
#' 
#' @return a function which takes in names and generates a plot with track views
#' and a summary plot per locus
#'
#' @export
#' @author Alicia Schep and Justin Finkle
#' 
#' @examples 
#' 
#' library(GenomicRanges)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' 
#' ## First we'll read in some sample chip-seq data
#' genomation_dir <- system.file("extdata", package = "genomationData")
#' samp.file <- file.path(genomation_dir,'SamplesInfo.txt')
#' samp.info <- read.table(samp.file, header=TRUE, sep='\t', 
#'                         stringsAsFactors = FALSE)
#' samp.info$fileName <- file.path(genomation_dir, samp.info$fileName)
#' 
#' ## we'll also read in some RNA counts
#' data(rpkm_chr21)
#' 
#' ## From the ranges of the rpkm object, we'll pull out the tss
#' tss <- promoters(SummarizedExperiment::rowRanges(rpkm_chr21),
#'                  up = 1, down = 1)
#' 
#' 
#' ## Make track plotter
#' 
#' track_plotter <- make_track_plotter(samp.info$fileName[1:3], 
#'   annotation = TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'   track_names = samp.info$sampleName[1:3] , 
#'   share_y = TRUE)
#'   
#' ## Make summary plotter
#' 
#' summary_plotter <- make_summary_plotter(rpkm_chr21,
#'   groups = "GROUP") 
#'   
#' ## Make combined plotter
#'   
#' combined_plotter <- make_track_plus_summary_plotter(track_plotter, 
#'   summary_plotter,
#'   resize(tss, width = 5000, fix = "center"),
#'   rownames(rpkm_chr21))
#'   
#' if (interactive()){
#'   combined_plotter(rownames(rpkm_chr21)[1:3])
#' }   
#' 
make_track_plus_summary_plotter <- function(track_function,
                                            summary_function,
                                            windows,
                                            row_names,
                                            summary_width = 0.25){
  out <- function(rows, track_args = NULL, summary_args = NULL){
    ix <- match(rows, row_names)
    tracks <- do.call(track_function, c(list(windows = windows[ix]), 
                                        track_args))
    summaries <- do.call(summary_function, c(list(row_names = row_names[ix]), 
                                             summary_args))
    new("GenomeTrackWidget", tracks = tracks, summaries = summaries, 
        summary_width = summary_width)
  }
  
  return(out)
}

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
              range <- c(min(object), max(object))
            } else{
              range <- NULL
            }
            
            if (length(object@xtitle) == 0){
              if (length(object) > 1 || object[[1]]@view@relative){
                xtitle <- "Relative Position"
              } else{
                xtitle <- as.character(seqnames(object[[1]]@view@range))
              }
            } else{
              xtitle <- object@xtitle
            }
            
            layout_setting <- 
              list(xaxis = 
                     list(title = xtitle,
                          zeroline = FALSE,
                          anchor = gsub("yaxis","y",
                                        ynames_flat[length(ynames_flat)]),
                          range = get_range(object[[1]]@view),
                          domain = x_domain))
            
            sizes <- unlist(purrr::map(as.list(object), function(y) y@heights ))
            
            sizes <- sizes / sum(sizes)
            
            domains <- list()
            start_domain <- 0
            k <- length(sizes)
            for (i in rev(seq_along(object))){
              domains[[i]] <-  list()
              for (j in rev(seq_along(object[[i]]))){
                domains[[i]][[j]] <- c(start_domain, start_domain +
                                         (sizes[k]*0.95))
                start_domain <- start_domain + sizes[k]
                k <- k - 1
              }
            }
            
            layout_setting <- c(layout_setting, 
                                unlist(purrr::pmap(list(object = 
                                                          as.list(object), 
                                                        yname = ynames,
                                                        domain = domains),
                                                   get_layout,
                                                   range = range), 
                                       recursive = FALSE))
            
            layout_setting
          })

setMethod(get_layout, signature = c(object = "LocusSummaryList"),
          definition = function(object, ynames, xax = "xaxis2", 
                                x_domain = c(0,1), ...){
            

            layout_setting <- list()
            layout_setting[[xax]] <- list(zeroline = FALSE,
                                          #showline = FALSE,
                                          anchor = gsub("yaxis","y",
                                                        ynames[length(ynames)]),
                                         domain = x_domain)
            
            sizes <- rep(1 / length(object), length(object))
            domains <- list()
            start_domain <- 0
            for (i in rev(seq_along(object))){
              domains[[i]] <- c(start_domain, start_domain + (sizes[i]*0.95))
              start_domain <- start_domain + sizes[i]
            }
            
            layout_setting <- c(layout_setting, 
                                unlist(purrr::pmap(list(object = 
                                                          as.list(object), 
                                                        yname = ynames,
                                                        domain = domains),
                                                   get_layout,
                                                   anchor = xax), 
                                                   recursive = FALSE))
            
            layout_setting
          })

get_range <- function(view){
  out <- c(relative_position(view, 
                      start(view@range)),
    relative_position(view, 
                      end(view@range)))
  if (out[2] < out[1]) out <- rev(out)
  out
}


multi_locus_to_plotly_list <- function(x){
  
  if (length(x@tracks) >= 1){
    
    lengths <- vapply(x@tracks, length, 0)
    track_ynames <- purrr::map2(as.list(x@tracks),
                                cumsum(lengths) - lengths[1] + 1,
                                yaxis_names)
    
    traces <- make_trace(x@tracks, track_ynames)
    
    if (length(x@summaries) == 0){
      x_domain <- c(0,1)
    } else{
      x_domain <- c(0, (1 - x@summary_width) * 0.95)
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
    layout_setting <- c(layout_setting, 
                        get_layout(x@summaries,
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
#' @rdname to_widget
setMethod(to_widget,
          c("LocusViewList"),
          function(p){
            p <- new("GenomeTrackWidget", tracks = p)
            out <- multi_locus_to_plotly_list(p)
            htmlwidgets::createWidget(
              name = "GenomicWidgets",
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
#' @rdname to_widget
setMethod(to_widget,
          c("LocusSummaryList"),
          function(p){
            p <- new("GenomeTrackWidget", summaries = p, summary_width = 1)
            out <- multi_locus_to_plotly_list(p)
            htmlwidgets::createWidget(
              name = "GenomicWidgets",
              x = out,
              width = out$layout$width,
              height = out$layout$height,
              sizingPolicy = htmlwidgets::sizingPolicy(browser.fill = TRUE,
                                                       viewer.fill = TRUE,
                                                       defaultWidth = "100%",
                                                       defaultHeight = 400),
              dependencies = plotly_dependency())
          })



#' to_widget
#' 
#' Method to convert GenomeTrackWidget to htmlwidgets objects
#' @param p GenomeTrackWidget or other object storing plot information
#' 
#' @return htmlwidgets object
#' @name to_widget
#' @rdname to_widget
#' @aliases to_widget,GenomeTrackWidget-method to_widget,NULL-method 
#' to_widget,LocusViewList-method to_widget,LocusView-method
#' @export
setMethod(to_widget,
          c("GenomeTrackWidget"),
          function(p){
            out <- multi_locus_to_plotly_list(p)
            htmlwidgets::createWidget(
              name = "GenomicWidgets",
              x = out,
              width = out$layout$width,
              height = out$layout$height,
              sizingPolicy = htmlwidgets::sizingPolicy(browser.fill = TRUE,
                                                       viewer.fill = TRUE,
                                                       defaultWidth = "100%",
                                                       defaultHeight = 400),
              dependencies = plotly_dependency())
          })

#' @rdname to_widget
setMethod(to_widget,
          c("NULL"),
          function(p){
            NULL
          })
