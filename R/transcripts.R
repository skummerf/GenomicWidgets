#' @export
#' @importClassesFrom GenomicFeatures TxDb
setMethod("unpack_transcripts", signature = "TxDb",
          function(object){
            unpack_transcripts_inner(object)
          })

#' @importClassesFrom OrganismDbi OrganismDb
setMethod("unpack_transcripts", signature = "OrganismDb",
          function(object){
            unpack_transcripts_inner(object)
          })

setMethod("unpack_transcripts", signature = "NULL",
          function(object){
            NULL
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
              tx_info <- add_stepping(tx_info, stacking = stacking)
            
            # Make AnnotationPlot object
            new("AnnotationPlot",
                transcripts = tx_info,
                trackname = name
            )
          })




#' @importMethodsFrom iheatmapr make_trace
#' @export
setMethod(make_trace, signature = c(x = "AnnotationPlot"),
          definition = function(x, yax, view, xax = "xaxis", ...){
            #Invisible Points
            anno_data <- as.data.frame(x@transcripts, row.names = NULL)
            if (nrow(anno_data) == 0) return(NULL)
            anno_data <- dplyr::mutate(anno_data,
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
            anno_data <- dplyr::group_by(anno_data, transcript)
            trace_data <- dplyr::summarise(anno_data, trace = list(list(x = midpoint,
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
            tx_info <- dplyr::mutate(tx_info,
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
            htmlwidgets::createWidget(name = "GenomicWidgets",
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

#' make_plotly_color
#' Convert R colors to plotly compatible rgb values
#'
#' @param color_str 
#'
#' @return
#' @export
#'
#' @author Justin Finkle
#' @examples
make_plotly_color <- function(color_str){
  # Plotly accepts some string colors, but it is safer to use rgb(0,0,0) values
  p_color <- rgb(t(col2rgb(color_str)), maxColorValue = 255)
  return(p_color)
}


#' make_rect
#' Make the rectangle objects. These are used to display annotation features such
#' as cds, and utr.
#'
#' @param df data.frame: data used to draw the shapes. Required columns include 
#' "start", "end", and "stepping"
#' @param height numeric: the height of the rectangles
#' @param yref character: axis on which to draw the rectangles, in the form of 
#' "yaxis[2-9]" or "y[2-9]"
#' @param fillcolor character: color of the rectangle
#'
#' @return list of rectangle shapes
#' @keywords internal
#'
#' @author Justin Finkle
make_rect <- function(df, 
                      height, 
                      yref,
                      fillcolor = 'lightslateblue'){
  fillcolor <- make_plotly_color(fillcolor)
  yref <- gsub("yaxis", "y", yref)
  if(nrow(df)>0){
    rect_list <- vector("list", nrow(df))
    for(e in 1:nrow(df)){
      row <- df[e, ]
      rect_list[[e]] <- list(type = "rect", fillcolor = fillcolor, opacity = 1, 
                             line=list(width=0),
                             x0 = row$start, x1 = row$end, xref = "x",
                             y0 = row$stepping-height, y1 = row$stepping+height, 
                             yref = yref)
    }
  } else {
    rect_list <- NULL
  }
  return(rect_list)
}

#' make_arrows
#' Draw arrows for introns
#'
#' @param df data.frame: data used to draw the shapes. Required columns include 
#' "start", "end", "strand", "midpoint", and "stepping"
#' @param yref character: axis on which to draw the rectangles, in the form of 
#' "yaxis[2-9]" or "y[2-9]"
#' @param arrowlen numeric: how long the arrow should be, in x coordinates
#' @param arrowheight numeric: how tall the arrow should be, in yref coordinates
#' @param arrowgap  numeric: gap between arrows on long introns
#'
#' @return list of arrow shapes
#'
#' @keywords interanl
#' @author Justin Finkle
make_arrows <- function(df, 
                        yref, 
                        arrowlen = 500, 
                        arrowheight = 0.15, 
                        arrowgap = 1500){
  # Set the y ref
  yref <- gsub("yaxis", "y", yref)
  if(nrow(df)>0){
    line_list <- vector("list", nrow(df))
    arrow_list <- vector("list", nrow(df))
    for(i in 1:nrow(df)){
      row <- df[i, ]
      
      # Make the intron line
      line_list[[i]] <- list(x0 = row$start, y0=row$stepping, 
                             x1 = row$end, y1 = row$stepping, 
                             xref = "x", yref = yref,
                             type = "line",
                             line = list(width = 0.5, 
                                         dash = ifelse(row$strand == "-", 
                                                       "dot","solid")))
      # Add arrows to the lines
      if (row$end - row$start > 2 * arrowlen){
        if (row$strand == "-"){
          arrow_pos <- row$midpoint - arrowlen * 0.5
        } else{
          arrow_pos <- row$midpoint + arrowlen * 0.5
        }
        arrow_list[[i]] <- arrow_helper(arrow_pos,
                                        row$strand,
                                        arrowlen,
                                        arrowheight,
                                        row$stepping,
                                        yref)
      }
    }
    out <- c(unlist(arrow_list, recursive = FALSE), line_list)
  } else {
    out <- NULL
  }
  return(out)
}


arrow_helper <- function(arrow_start, 
                         strand, 
                         arrowlen, 
                         arrowheight, 
                         y, 
                         yref){
  arrow_end <- ifelse(strand =="-", arrow_start + arrowlen, 
                      arrow_start - arrowlen)
  list(list(x0 = arrow_start, x1=arrow_end, 
            y0 = y, y1 = y - arrowheight, 
            xref = "x", yref = yref,
            type = "line",
            line = list(width = 0.5)),
       list(x0 = arrow_start, x1=arrow_end, 
            y0 = y, y1 = y + arrowheight, 
            xref = "x", yref = yref,
            type = "line",
            line = list(width = 0.5)))
}



