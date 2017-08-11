setMethod("unpack_transcripts", signature = "TxDb",
          function(object){
            unpack_transcripts_inner(object)
          })



unpack_transcripts_inner <- function(db_object){
  message("Loading transcript annotation data...")
  tx_data <- new("TranscriptParts",
                 intron = intronsByTranscript(db_object, use.names=TRUE),
                  utr5 = fiveUTRsByTranscript(db_object, use.names=TRUE),
                  utr3 = threeUTRsByTranscript(db_object, use.names=TRUE),
                  cds = cdsBy(db_object, by='tx', use.names=TRUE),
                  exon = exonsBy(db_object, by='tx', use.names=TRUE),
                  transcripts = transcripts(db_object, columns = c("TXID","TXNAME")),
                  seqlevelsStyle =  seqlevelsStyle(db_object))
  return(tx_data)
}

#TO DO: ADD PRINT METHOD FOR TRANSCRIPTS

setMethod("subset_transcripts", signature = c("GRanges","TxDb"),
          function(window, object, no_introns = FALSE, ...){
            parts <- unpack_transcripts(object)
            subset_transcripts(window, parts, no_introns = no_introns, ...)
          })


setMethod("subset_transcripts", signature = c("GRanges","TranscriptParts"),
          function(window, object, no_introns = FALSE, ...){
            in_style <- seqlevelsStyle(window)[[1]]
            seqlevelsStyle(window) <- object@seqlevelsStyle[1]
            tx <- subsetByOverlaps(object@transcripts, window, maxgap = 0L,
                                   minoverlap = 1L, type = "any")
            tx_names <- unlist(tx$TXNAME)
            gr <- get_tx_parts(tx_names, object)
            if(length(gr)){
              seqlevelsStyle(gr) <- in_style
              if(no_introns){
                gr <- gr[mcols(gr)[['feature']]!='intron']
              }
            }
            # Some tx only partially overlap and some tx parts outside the window get included
            gr <- subsetByOverlaps(gr, window)
            return(gr)
          })


get_tx_parts <- function(tx_names, tx_data){
  if(!length(tx_names)){ return(GRanges())}
  parts <- GenomicRanges::GRanges()
  for (n in c("intron","utr5","utr3","cds","exon")){
    tx_subset <- tx_names[tx_names %in% names(slot(tx_data,n))]
    if (length(tx_subset)){
      part_gr <- unlist(slot(tx_data,n)[tx_subset])
      if (length(part_gr)){
        part_gr$transcript <- names(part_gr)
        part_gr$feature <- n
        if(!('exon_name' %in% colnames(mcols(part_gr)))){
          part_gr$exon_name <- NA
        }
        part_gr <- part_gr[, c('transcript', 'feature', 'exon_name')]
        parts <- c(parts, part_gr)
      }
    }
  }

  # Find ncRNA in exons
  exon_tx <- unique(parts$transcript[parts$feature=='exon'])
  coding_tx <- unique(parts$transcript[parts$feature!='intron' & parts$feature!='exon'])

  # ncRNA are exons that aren't in CDS or UTRs
  nc_tx <- setdiff(exon_tx, coding_tx)
  parts$feature[parts$transcript %in% nc_tx & parts$feature=='exon'] <-'ncRNA'

  # Remove the now redudant exons
  parts <- parts[parts$feature != 'exon']

  return(parts)
}


setMethod("add_stepping", signature = c("GRanges"),
          function(object, stacking = c('squish', 'dense'), ...){
            stacking <- match.arg(stacking)
            if(length(object)){
              if(stacking == 'squish'){
                stopifnot("transcript" %in% colnames(mcols(object)))
                obj.lst <- split(object, as.character(seqnames(object)))
                lv <- endoapply(obj.lst, function(x) {
                      x.n <- split(x, values(x)[, "transcript"])
                      irs <- unlist(range(ranges(x.n)))
                      irs.new <- resize(irs, fix = "center", width = width(irs))
                      irs.new <- sort(irs.new)
                      .lvs <- disjointBins(irs.new)
                      values(x)$stepping <- .lvs[as.character(values(x)[,"transcript"])]
                      x
                })
                object <- unlist(lv)
              } else if(stacking == 'dense'){
                # The simple solution for now
                mcols(object)$stepping <- 1
              }}
              return(object)
            })

setMethod("make_annotation_track", c("GRanges", "TxDb"),
          function(window, object, stacking = c('squish', 'dense'), ...){
            
            make_annotation_track(as_view_range(window), 
                                  unpack_transcripts(object), 
                                  match.arg(stacking), ...)
            
          })

setMethod("make_annotation_track", c("ViewRange", "TxDb"),
          function(window, object, stacking = c('squish', 'dense'), ...){
            
            make_annotation_track(window, 
                                  unpack_transcripts(object), 
                                  match.arg(stacking), ...)
            
          })

setMethod("make_annotation_track", c("GRanges", "TranscriptParts"),
          function(window, object, stacking = c('squish', 'dense'), ...){

          make_annotation_track(as_view_range(window), object, 
                                match.arg(stacking), ...)
          
          })

setMethod("make_annotation_track", c("ViewRange", "TranscriptParts"),
          function(window, object, stacking = c('squish', 'dense'), 
                   name = "", ...){
            
            stacking <- match.arg(stacking)
            
            # subset transcript info
            tx_info <- subset_transcripts(window, object)
            
            # add stepping
            if (length(tx_info) >= 1)
              tx_info <- add_tx_stepping(tx_info, stacking = stacking)
            
            # Make AnnotationPlot object
            new("AnnotationPlot",
                transcripts = tx_info,
                trackname = name
            )
          })




#' @importMethodsFrom iheatmapr make_trace
#' @import dplyr
#' @export
setMethod(make_trace, signature = c(x = "AnnotationPlot"),
          definition = function(x, yax, view, xax = "xaxis", ...){
            #Invisible Points
            anno_data <- as.data.frame(x@transcripts, row.names = NULL)
            if (nrow(anno_data) == 0) return(NULL)
            anno_data <- mutate(anno_data,
                               text = paste0('Tx ID: ',
                                     transcript,
                                     '<br>',
                                     feature,
                                     '<br>',
                                     "strand: ",
                                     strand),
                               start = relative_position(view, start),
                               end = relative_position(view, end),
                               midpoint = (start + end) / 2)
            anno_data <- group_by(anno_data, transcript)
            trace_data <- summarise(anno_data, trace = list(list(x = midpoint,
                                                  y = stepping,
                                                  text = text,
                                                  yaxis = gsub("yaxis","y",yax),
                                                  xaxis = gsub("xaxis","x",xax),
                                                  hoverinfo = 'x+text',
                                                  opacity = 0,
                                                  type='scatter',
                                                  showlegend = FALSE,
                                                  name = unique(transcript),
                                                  mode = 'markers')))
            trace_data$trace
          })

setMethod("make_shapes", c(x = "AnnotationPlot"),
          function(x, yax, view, ...){
            ann_ax <- gsub("yaxis","y",yax)
            tx_info <- as.data.frame(x@transcripts, row.names = NULL)
            tx_info <- mutate(tx_info,
                              start = relative_position(view, start),
                              end = relative_position(view, end),
                              midpoint = (start + end) / 2)

            cds_rect <- make_rect(tx_info[tx_info$feature == 'cds', ],
                                  height = 0.4,
                                  ann_ax)
            utr_rect <- make_rect(tx_info[grep("utr", tx_info$feature), ],
                                  height=0.25,
                                  ann_ax)
            ncRNA_rect <- make_rect(tx_info[tx_info$feature == 'ncRNA', ],
                                    height=0.25,
                                    ann_ax)
            intron_arrow <- make_arrows(tx_info[tx_info$feature == 'intron', ], ann_ax,
                                        arrowlen = width(view) * 0.01)
            # Compile new shapes
            tx_shapes <- c(cds_rect, utr_rect, ncRNA_rect, intron_arrow)

            return(tx_shapes)
          })




annotation_to_plotly_list <- function(x){
  traces <- make_trace(x, "yaxis")
  shapes <- make_shapes(x, "yaxis")
  layout_setting <- list()#get_layout(x)
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
          c("AnnotationPlot"),
          function(p){
            out <- annotation_to_plotly_list(p)
            print("a")
            htmlwidgets::createWidget(name = "chipVis",
                         x = out,
                         width = out$layout$width,
                         height = out$layout$height,
                         sizingPolicy = htmlwidgets::sizingPolicy(browser.fill = TRUE,
                                                                  viewer.fill = TRUE,
                                                     defaultWidth = "100%",
                                                     defaultHeight = 400),
                        dependencies = plotly_dependency())
          })



plotly_dependency <- function(){
  htmltools::htmlDependency(
    "plotlyjs", "1.29.2",
    src = system.file('htmlwidgets', 'lib', 'plotlyjs', package = 'iheatmapr'),
    script = "plotly-latest.min.js",
    stylesheet = "plotly-htmlwidgets.css"
  )
}




