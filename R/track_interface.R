#' set_track_parameters
#' 
#' Setup parameters for plotting coverage signals along genome tracks.
#' Result from this function can be passed to \code{\link{plot_tracks}}
#'
#' @param object vector of bam or bigwig file names
#' @param annotation TxDb or OrganismDb object
#' @param ... additional arguments
#' @param track_names names to associate with each file
#' @param groups vector of group assignments.  traces will be grouped onto 
#' subplots based on group assignments (if only showing 1 region)
#' @param share_y share the y axis?
#' @param fill fillmode for line plot
#' @param showlegend show the legend?
#' @param colors colors for each bam file
#' @param mode mode for plot
#' @param annotation_position plot annotations on bottom or on top of signal 
#' traces
#' @param annotation_size relative size of annotation plot
#' @param summary Summary parameters from \code{\link{set_summary_parameters}}
#' @param layout list of additional plotly layout arguments
#' @export
#' @rdname set_track_parameters
#' @name set_track_parameters
#' @aliases set_track_parameters,character-method
#' 
#' @author Alicia Schep and Justin Finkle
#' @return object storing track parameters, for use in \code{\link{plot_tracks}}
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
#' track_params <- set_track_parameters(samp.info$fileName[1:3], 
#'   annotation = TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'   track_names = samp.info$sampleName[1:3] , 
#'   share_y = TRUE)
#'   
#' if (interactive()){
#'   plot_tracks(ctcf.peaks[1], track_params)
#'   plot(ctcf.peaks[1:3], track_params)
#' }   
#' 
setMethod("set_track_parameters", c("character"),
          function(object,
                   annotation = NULL,
                   track_names = if (!is.null(names(object)))
                     names(object) 
                   else basename(object),
                   groups = track_names,
                   share_y = TRUE,
                   showlegend = TRUE,
                   colors = NULL,
                   fill = c("tozeroy","none"),
                   mode = "lines",
                   annotation_position = c("bottom","top"),
                   annotation_size = 0.25,
                   summary = NULL,
                   layout = list()){
            
            params <- list(Class = "TrackParameters",
                           data = object,
                           annotation = unpack_transcripts(annotation),
                           track_names = track_names,
                           groups = as.factor(groups),
                           summary = summary,
                           showlegend = showlegend,
                           share_y = share_y,
                           colors = select_colors(colors,length(object)),
                           fill = match.arg(fill),
                           annotation_position = match.arg(annotation_position),
                           mode = mode,
                           annotation_size = annotation_size,
                           layout = layout)
            
            do.call(new, params)
            
          })

#' set_summary_parameters
#' 
#' Setup parameters for plotting summaries along genome tracks.
#' Result from this function can be passed to \code{\link{set_track_parameters}}
#' 
#' @param object SummarizedExperiment
#' @param assay_name name of assay to use
#' @param ... additional arguments
#' @param groups either vector of group assignments of name of column in object 
#' colData that corresponds to vector of group assignments
#' @param showlegend show the legend?
#' @param colors colors to use
#' @param boxpoints plot individual points?
#' @param pointpos relative position of points to boxes
#' @param ytitle name for yaxis
#' @param ranges ranges corresponding to rows of object
#' @param width relative width of summary plots when plotting tracks
#' @export
#' @return object storing summary parameters, for use in
#'  \code{\link{set_track_parameters}}
#' @author Alicia Schep and Justin Finkle
#' @rdname set_summary_parameters
#' @name set_summary_parameters
#' @aliases set_summary_parameters,SummarizedExperiment-method
#' set_summary_parameters,RangedSummarizedExperiment-method
#' @examples 
#' 
#' 
#' library(GenomicRanges)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' 
#' ## we'll read in some RNA counts
#' data(rpkm_chr21)
#' 
#' ## From the ranges of the rpkm object, we'll pull out the tss
#' chr21_promoters <- promoters(SummarizedExperiment::rowRanges(rpkm_chr21),
#'                  up = 1000, down = 1000)
#' 
#' ## set summary parameters
#' 
#' summary_params <- set_summary_parameters(rpkm_chr21,
#'   groups = "GROUP", ranges = chr21_promoters) 
#'   
#' ## We'll also read in some track data to plot
#' genomation_dir <- system.file("extdata", package = "genomationData")
#' samp.file <- file.path(genomation_dir,'SamplesInfo.txt')
#' samp.info <- read.table(samp.file, header=TRUE, sep="\t", 
#'                         stringsAsFactors = FALSE)
#' samp.info$fileName <- file.path(genomation_dir, samp.info$fileName)
#' 
#' ## Make track plotter using summary parametrs
#' 
#' track_params <- set_track_parameters(samp.info$fileName[1:3], 
#'   annotation = TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'   track_names = samp.info$sampleName[1:3], 
#'   share_y = TRUE,
#'   summary = summary_params)
#'   
#' if (interactive()){
#'   plot_tracks(rownames(rpkm_chr21)[1:3], track_params)
#'   plot(chr21_promoters[1:3], track_params)
#' }   
#' 
setMethod("set_summary_parameters", c("SummarizedExperiment"),
          function(object,  
                   ranges,
                   assay_name = assayNames(object)[1],
                   groups = colnames(object),
                   colors = "blue",
                   showlegend = length(colors) > 1, 
                   boxpoints = c("all","Outliers","false"),
                   pointpos = 0,
                   ytitle = "Expression",
                   width = 0.3){
            
            # Check groups
            if (is.null(groups)){
              groups <- ""
            } else if (length(groups) == 1){
              if (groups %in% colnames(colData(object))) 
                groups <- colData(object)[,groups]
            } else{
              stopifnot(length(groups) == ncol(object))
            }
            
            params <- list(Class = "SummaryParameters",
                           data = object,
                           assay_name = assay_name,
                           groups = groups,
                           colors = colors,
                           showlegend = showlegend,
                           boxpoints = match.arg(boxpoints),
                           pointpos = pointpos,
                           ytitle = ytitle,
                           width = width,
                           ranges = ranges
                           )
            
            do.call(new, params)
            
          })

#' @export
#' @rdname set_summary_parameters
setMethod("set_summary_parameters", c("RangedSummarizedExperiment"),
          function(object, 
                   ranges = rowRanges(object),
                   assay_name = assayNames(object)[1],
                   groups = NULL,
                   colors = "blue",
                   showlegend = length(colors) > 1, 
                   boxpoints = c("all","Outliers","false"),
                   pointpos = 0,
                   ytitle = "Expression",
                   width = 0.3){
            
            # Check groups
            if (is.null(groups)){
              groups <- ""
            } else if (length(groups) == 1){
              if (groups %in% colnames(colData(object))) 
                groups <- colData(object)[,groups]
            } else{
              stopifnot(length(groups) == ncol(object))
            }
            
            params <- list(Class = "SummaryParameters",
                           data = object,
                           assay_name = assay_name,
                           groups = groups,
                           colors = colors,
                           showlegend = showlegend,
                           boxpoints = match.arg(boxpoints),
                           pointpos = pointpos,
                           ytitle = ytitle,
                           width = width,
                           ranges = ranges
            )
            
            do.call(new, params)
          })

#' @export
#' @rdname plot_tracks
setMethod("plot_tracks", c("GenomicRanges"),
          function(windows, 
                   params, 
                   locus_names = mcols(windows)$name,
                   offset = width(windows[1]) %/% 2, 
                   xtitle = 
                     if (length(windows) > 1) "Relative Position" 
                      else seqnames(windows),
                   ..., 
                   summary_args = list()){
            
            
            arglist <- modify_param_args(params, xtitle, offset, locus_names,
                                         ...)
            tracks <- do.call(multi_locus_view, 
                              c(list(windows = windows), arglist))
            
            out <- new("GenomeTrackWidget", 
                       tracks = tracks, 
                       summary_width = 0,
                       layout = arglist$layout)
            
            if (!is.null(params@summary)){
              ix <- match(windows, params@summary@ranges)
              if (any(is.na(ix))){
                warning("Windows don't match ranges stored in TrackParameters ",
                        "for summary.  Not plotting summaries.")
              } else{
                
                summary_arglist <- modify_summary_args(params,summary_args)
                
                summary_arglist[["row_names"]] <- 
                  rownames(params@summary@data)[ix]
                
                summaries <- do.call(make_locus_summaries, summary_arglist)
                out <- new("GenomeTrackWidget",
                           tracks = tracks, 
                           summaries = summaries, 
                           summary_width = params@summary@width,
                           layout = arglist$layout)
              }
            } else{
              out <- new("GenomeTrackWidget",
                         tracks = tracks, 
                         summary_width = 0,
                         layout = arglist$layout)
            }
            out
          })

#' plot_tracks
#' 
#' @param windows GenomicRanges or rownames for summary object
#' @param params TrackParameters, from \code{\link{set_track_parameters}}
#' @param locus_names names for each genomic locus represented by windows
#' @param offset offset to use for center of region, used when plotting multiple 
#' regions
#' @param xtitle title for x axis
#' @param ... additional arguments from \code{\link{set_track_parameters}}
#' @param summary_args lof arguments to override  from 
#' \code{\link{set_summary_parameters}}
#' @export
#' @aliases plot_tracks,character-method plot_tracks,GenomicRanges-method
#' @author Alicia Schep and Justin Finkle
#' @rdname plot_tracks
#' @name plot_tracks
#' @return GenomeTrackWidgets object, which is displayed as htmlwidgets. To 
#' convert manually to htmlwidgets, use to_widget
#' 
#' @examples 
#' 
#' library(GenomicRanges)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' 
#' ## we'll read in some RNA counts
#' data(rpkm_chr21)
#' 
#' ## From the ranges of the rpkm object, we'll pull out the tss
#' chr21_promoters <- promoters(SummarizedExperiment::rowRanges(rpkm_chr21),
#'                  up = 1000, down = 1000)
#' 
#' ## set summary parameters
#' 
#' summary_params <- set_summary_parameters(rpkm_chr21,
#'   groups = "GROUP", ranges = chr21_promoters) 
#'   
#' ## We'll also read in some track data to plot
#' genomation_dir <- system.file("extdata", package = "genomationData")
#' samp.file <- file.path(genomation_dir,'SamplesInfo.txt')
#' samp.info <- read.table(samp.file, header=TRUE, sep="\t", 
#'                         stringsAsFactors = FALSE)
#' samp.info$fileName <- file.path(genomation_dir, samp.info$fileName)
#' 
#' ## Make track plotter using summary parametrs
#' 
#' track_params <- set_track_parameters(samp.info$fileName[1:3], 
#'   annotation = TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'   track_names = samp.info$sampleName[1:3], 
#'   share_y = TRUE,
#'   summary = summary_params)
#'   
#' if (interactive()){
#'   plot_tracks(rownames(rpkm_chr21)[1:3], track_params)
#'   plot(chr21_promoters[1:3], track_params)
#' }   
#' 
setMethod("plot_tracks", c("character"),
          function(windows, 
                   params, 
                   locus_names = windows,
                   offset = NULL, 
                   xtitle = 
                     if (length(windows) > 1) "Relative Position" 
                      else NULL,
                   ..., 
                   summary_args = list()){
            
            
            
            if (is.null(params@summary)){
              stop("params must include summary parameters for calling ", 
                   "plot_tracks via a character vector")
            }
            
            
            ix <- match(windows, rownames(params@summary@data))
            ranges <- params@summary@ranges[ix]
            
            if (is.null(offset)){offset <- width(ranges[1]) %/% 2}
            if (is.null(xtitle)){
              if (length(windows) > 1){
                xtitle <- ""
              } else{
                xtitle <- seqnames(ranges)
              }
            }
            
            arglist <- modify_param_args(params, xtitle, offset, locus_names,
                                         ...)
            tracks <- do.call(multi_locus_view, 
                              c(list(windows = ranges), arglist))
            
            
            if (any(is.na(ix))){
              warning("Windows don't match ranges stored in TrackParameters ","
                        for summary.  Not plotting summaries.")
            } else{
              
              summary_arglist <- modify_summary_args(params,summary_args)
              
              summary_arglist[["row_names"]] <- 
                windows
              
              summaries <- do.call(make_locus_summaries, summary_arglist)
              out <- new("GenomeTrackWidget", tracks = tracks, 
                         summaries = summaries, 
                         summary_width = params@summary@width,
                         layout = arglist$layout)
            }
            out
          })

modify_param_args <- function(params, xtitle, offset, locus_names, ...){
  
  new_args <- list(...)
  
  new_layout <- if ("layout" %in% names(new_args)){
    modifyList(params@layout, new_args[["layout"]])
  } else{
    params@layout
  }
  
  default_arglist <- list(
    object = params@data,
    annotation = params@annotation,
    track_names = params@track_names,
    name = locus_names,
    groups = params@groups,
    share_y = params@share_y,
    fill = params@fill,
    showlegend = params@showlegend,
    colors = params@colors,
    mode = params@mode,
    annotation_position = params@annotation_position,
    annotation_size = params@annotation_size,
    layout = new_layout,
    offset = offset,
    xtitle = xtitle
  )
  
  modifyList(default_arglist, new_args)
  
}

modify_summary_args <- function(params, summary_args){
  
  default_summary_arglist <- 
    list(object = params@summary@data,
         assay_name = params@summary@assay_name,
         groups = params@summary@groups,
         showlegend = params@summary@showlegend,
         boxpoints = params@summary@boxpoints,
         pointpos = params@summary@pointpos,
         ytitle = params@summary@ytitle)
  
  modifyList(default_summary_arglist, 
                                summary_args)
  
}
