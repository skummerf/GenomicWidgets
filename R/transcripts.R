#' unpack_transcripts
#' 
#' Method to load transcript annotation data to a format that is easy to make
#' transcript visualizations from. The resulting object can be passed to 
#' make_track_plotter and will result in faster track function creation than 
#' passing the TxDb or OrganismDb object. Useful if calling make_track_plotter
#' multiple times.
#' @param object a TxDb or OrganismDb object
#' @export
#' @name unpack_transcripts
#' @rdname unpack_transcripts
#' @return TranscriptParts object
#' @author Alicia Schep and Justin Finkle
#' @importClassesFrom GenomicFeatures TxDb
#' @import GenomicFeatures
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @aliases unpack_transcripts,TxDb-method 
#' unpack_transcripts,OrganismDbi-method unpack_transcripts,NULL-method
#' unpack_transcripts,TranscriptParts-method
#' @examples 
#' 
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' transc <- unpack_transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' class(transc)
setMethod("unpack_transcripts", signature = "TxDb",
          function(object){
            unpack_transcripts_inner(object)
          })

#' @export
#' @rdname unpack_transcripts
#' @importClassesFrom OrganismDbi OrganismDb
setMethod("unpack_transcripts", signature = "OrganismDb",
          function(object){
            unpack_transcripts_inner(object)
          })

setMethod("unpack_transcripts", signature = "NULL",
          function(object){
            NULL
          })

setMethod("unpack_transcripts", signature = "TranscriptParts",
          function(object){
            object
          })


unpack_transcripts_inner <- function(db_object){
  message("Loading transcript annotation data...")
  tx_data <- new("TranscriptParts",
                 intron = intronsByTranscript(db_object, use.names=TRUE),
                  utr5 = fiveUTRsByTranscript(db_object, use.names=TRUE),
                  utr3 = threeUTRsByTranscript(db_object, use.names=TRUE),
                  cds = cdsBy(db_object, by="tx", use.names=TRUE),
                  exon = exonsBy(db_object, by="tx", use.names=TRUE),
                  transcripts = transcripts(db_object, 
                                            columns = c("TXID","TXNAME")),
                  seqlevelsStyle =  seqlevelsStyle(db_object))
  return(tx_data)
}

#TO DO: ADD PRINT METHOD FOR TRANSCRIPTS

setMethod("subset_transcripts", signature = c("GRanges","TxDb"),
          function(window, object, no_introns = FALSE, ...){
            parts <- unpack_transcripts(object)
            subset_transcripts(window, parts, no_introns = no_introns, ...)
          })

#' @importFrom IRanges subsetByOverlaps IRanges
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
                gr <- gr[mcols(gr)[["feature"]]!="intron"]
              }
            }
            # some tx parts outside the window get included
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
        if(!("exon_name" %in% colnames(mcols(part_gr)))){
          part_gr$exon_name <- NA
        }
        part_gr <- part_gr[, c("transcript", "feature", "exon_name")]
        parts <- c(parts, part_gr)
      }
    }
  }

  # Find ncRNA in exons
  exon_tx <- unique(parts$transcript[parts$feature=="exon"])
  coding_tx <- unique(parts$transcript[parts$feature!="intron" & 
                                         parts$feature!="exon"])

  # ncRNA are exons that aren't in CDS or UTRs
  nc_tx <- setdiff(exon_tx, coding_tx)
  parts$feature[parts$transcript %in% nc_tx & parts$feature=="exon"] <- "ncRNA"

  # Remove the now redudant exons
  parts <- parts[parts$feature != "exon"]

  return(parts)
}


setMethod("add_stepping", signature = c("GRanges"),
          function(object, stacking = c("squish", "dense"), ...){
            stacking <- match.arg(stacking)
            if(length(object)){
              if(stacking == "squish"){
                stopifnot("transcript" %in% colnames(mcols(object)))
                obj.lst <- split(object, as.character(seqnames(object)))
                lv <- endoapply(obj.lst, function(x) {
                      x.n <- split(x, values(x)[, "transcript"])
                      irs <- unlist(range(ranges(x.n)))
                      irs.new <- resize(irs, fix = "center", width = width(irs))
                      irs.new <- sort(irs.new)
                      .lvs <- disjointBins(irs.new)
                      values(x)$stepping <- 
                        .lvs[as.character(values(x)[,"transcript"])]
                      x
                })
                object <- unlist(lv)
              } else if(stacking == "dense"){
                # The simple solution for now
                mcols(object)$stepping <- 1
              }}
              return(object)
            })

setMethod("make_annotation_track", c("GRanges", "TxDb"),
          function(window, object, stacking = c("squish", "dense"), ...){
            
            make_annotation_track(window, 
                                  unpack_transcripts(object), 
                                  match.arg(stacking), ...)
            
          })



setMethod("make_annotation_track", c("GRanges", "TranscriptParts"),
          function(window, object, stacking = c("squish", "dense"), 
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
setMethod(make_trace, signature = c(x = "AnnotationPlot"),
          definition = function(x, yax, view, xax = "xaxis", ...){
            #Invisible Points
            anno_data <- as.data.frame(x@transcripts, row.names = NULL)
            if (nrow(anno_data) == 0) return(NULL)
            anno_data$text <- paste0("Tx ID: ",
                                     anno_data$transcript,
                                     "<br>",
                                     anno_data$feature,
                                     "<br>",
                                     "strand: ",
                                     anno_data$strand)
            anno_data$start <- relative_position(view, anno_data$start)
            anno_data$end <- relative_position(view, anno_data$end)
            anno_data$midpoint <- (anno_data$start + anno_data$end) / 2
            
            base_list <- list(yaxis = gsub("yaxis","y",yax),
                              xaxis = gsub("xaxis","x",xax),
                              hoverinfo = "x+text",
                              opacity = 0,
                              type = "scatter",
                              showlegend = FALSE,
                              mode = "markers")
            traces <- lapply(unique(anno_data$transcript), function(tname){
              ix <- which(anno_data$transcript == tname)
              c(base_list, list(x = anno_data$midpoint[ix],
                   y = anno_data$stepping[ix],
                   name = tname))
              
            })
            traces
          })

setMethod("make_shapes", c(x = "AnnotationPlot"),
          function(x, yax, view, ...){
            ann_ax <- gsub("yaxis","y",yax)
            tx_info <- as.data.frame(x@transcripts, row.names = NULL)
            tx_info$start <- relative_position(view, tx_info$start)
            tx_info$end <- relative_position(view, tx_info$end)
            tx_info$midpoint <- (tx_info$start + tx_info$end) / 2

            cds_rect <- make_rect(tx_info[tx_info$feature == "cds", ],
                                  height = 0.4,
                                  ann_ax)
            utr_rect <- make_rect(tx_info[grep("utr", tx_info$feature), ],
                                  height=0.25,
                                  ann_ax)
            ncRNA_rect <- make_rect(tx_info[tx_info$feature == "ncRNA", ],
                                    height=0.25,
                                    ann_ax)
            intron_arrow <- make_arrows(tx_info[tx_info$feature == "intron", ],
                                        ann_ax,
                                        arrowlen = width(view@range) * 0.01)
            # Compile new shapes
            tx_shapes <- c(cds_rect, utr_rect, ncRNA_rect, intron_arrow)

            return(tx_shapes)
          })




#' make_plotly_color
#' Convert R colors to plotly compatible rgb values
#'
#' @param color_str 
#'
#' @return color rgb value
#'
#' @keywords internal
#' @author Justin Finkle
make_plotly_color <- function(color_str){
  # Plotly accepts some string colors, but it is safer to use rgb(0,0,0) values
  p_color <- rgb(t(col2rgb(color_str)), maxColorValue = 255)
  return(p_color)
}


#' make_rect
#' 
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
                      fillcolor = "lightslateblue"){
  fillcolor <- make_plotly_color(fillcolor)
  yref <- gsub("yaxis", "y", yref)
  if(nrow(df)>0){
    rect_list <- vector("list", nrow(df))
    for(e in seq_len(nrow(df))){
      row <- df[e, ]
      rect_list[[e]] <- list(type = "rect", 
                             fillcolor = fillcolor, 
                             opacity = 1, 
                             line=list(width=0),
                             x0 = row$start, 
                             x1 = row$end, 
                             xref = "x",
                             y0 = row$stepping-height, 
                             y1 = row$stepping+height, 
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
    for(i in seq_len(nrow(df))){
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



