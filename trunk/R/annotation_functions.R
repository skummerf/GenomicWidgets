# =============================================================================
# =============================================================================
# Transcript annotation functions
#' Load Transcript Data
#'
#' @param db_object database object: currently accepts TxDB and OrganismDB
#'
#' @return
#' @export
#'
#' @examples
load_tx_data <- function(db_object){
  db_class <- get_db_class(db_object)
  message("Loading transcript annotation data...")
  tx_data <- list(intron = intronsByTranscript(db_object, use.names=TRUE),
                  utr5 = fiveUTRsByTranscript(db_object, use.names=TRUE),
                  utr3 = threeUTRsByTranscript(db_object, use.names=TRUE),
                  cds = cdsBy(db_object, by='tx', use.names=TRUE),
                  exon = exonsBy(db_object, by='tx', use.names=TRUE))
  return(tx_data)
}

transcriptsByOverlaps <- function(x, ranges, maxgap = 0L, minoverlap = 1L,
                                  type = c("any", "start", "end"),
                                  columns = c("TXID", "TXNAME")) {
  subsetByOverlaps(transcripts(x, columns = columns), ranges, maxgap = maxgap, 
                   minoverlap = minoverlap, type = match.arg(type))
}

#' Title
#'
#' @param db_object 
#' @param range 
#' @param tx_data 
#' @param no_introns 
#'
#' @return
#' @export
#'
#' @examples
get_tx_annotation <- function(db_object, range, tx_data, no_introns=FALSE){
  in_style <- seqlevelsStyle(range)[[1]]
  seqlevelsStyle(range) <- seqlevelsStyle(db_object)
  tx <- chipVis:::transcriptsByOverlaps(db_object, range)
  tx_names <- unlist(tx$TXNAME)
  gr <- get_tx_features(tx_names, tx_data)
  if(length(gr)){
    seqlevelsStyle(gr) <- in_style
    if(no_introns){
      gr <- gr[mcols(gr)[['feature']]!='intron']
    }
  }
  # Some tx only partially overlap and some tx parts outside the range get included
  gr <- subsetByOverlaps(gr, range)
  return(gr)
}

get_tx_features <- function(tx_names, tx_data){
  if(!length(tx_names)){ return(GRanges())}
  for(n in names(tx_data)){
    if(length(tx_data[[n]])){
      tx_subset <- tx_names[tx_names %in% names(tx_data[[n]])]
      part_gr <- unlist(tx_data[[n]][tx_subset])
      part_gr$transcript <- names(part_gr)
      part_gr$feature <- n
      if(!('exon_name' %in% colnames(mcols(part_gr)))){
        part_gr$exon_name <- NA
      }
      part_gr <- part_gr[, c('transcript', 'feature', 'exon_name')]
    } else {
      part_gr <- GRanges()
    }
    if(exists("parts")){
      parts <- c(parts, part_gr)
    } else {
      parts <- part_gr
    }
  }
  
  # Find ncRNA in exons
  exon_tx <- unique(parts$transcript[parts$feature=='exon'])
  coding_tx <- unique(parts$transcript[parts$feature!='intron' & parts$feature!='exon'])
  
  # ncRNA are exons that aren't in CDS or UTRs
  nc_tx <- setdiff(exon_tx, coding_tx)
  parts$feature[parts$transcript %in% nc_tx & parts$feature=='exon'] <-'ncRNA'
  
  return(parts)
}