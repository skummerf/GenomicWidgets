#' @export
setAs("ViewRange","GRanges",
      function(from){
        out <- GRanges(seqnames = seqnames(from),
                strand = strand(from),
                ranges = IRanges(start = start(from), end = end(from)),
                seqinfo = seqinfo(from))
        mcols(out) <- mcols(from)
        return(out)
      })

#' @export
setAs("GRanges","ViewRange",
      function(from){
        if (!("name" %in% colnames(mcols(from)))) 
          mcols(from)$name <- as.character(seq_along(from))
        new("ViewRange", from)
      })

#' @export
setAs("GRanges","RelativeViewRange",
      function(from){
        if (!("name" %in% colnames(mcols(from)))) 
          mcols(from)$name <- ""
        if (!("relative" %in% colnames(mcols(from)))) 
          mcols(from)$relative <- start(from) + width(from) %/% 2
        new("RelativeViewRange", from)
      })

#' make_view_range
#' 
#' Makes a "ViewRange" object
#' 
#' @import BiocGenerics
#' @importFrom IRanges IRanges
#' @export
make_view_range <- function(chrom,
                       start, 
                       end, 
                       strand = c("*", "+", "-"),
                       name = character(0),
                       relative = FALSE){
  
  strand <- match.arg(strand)
  gr <- GRanges(Rle(c(chrom)), IRanges(start = start, end = end), 
                strand = strand)
  if (!relative){
    out <- new("ViewRange", 
               gr, 
               name = ifelse(length(name) == 0, as.character(gr), name))
  } else{
    out <- new("RelativeViewRange", 
               gr, 
               name = ifelse(length(name) == 0, as.character(gr), name))
  }
  out
}

#' fetch_view_range
#' 
#' Function to get a range around a gene, using an OrganismDb or TxDb object
#' 
#' @param db_object OrganismDb or TxDb object
#' @param symbol name of gene
#' @param keytype type of key
#' @param relative get range relative to TSS, TTS, or full gene
#' @param up basepairs upstream to include
#' @param down basepairs downstream to include
#' 
#' @return GenomicRanges
#' @author Alicia Schep and Justin Finkle
#' @export
#' 
#' @examples 
#' 
#' library(Homo.sapiens)
#' fetch_view_range(Homo.sapiens,"GLI2")
fetch_view_range <- function(db_object, 
                 symbol, 
                 keytype = c('SYMBOL', 'TXNAME'),
                 relative = c("TSS", "TTS", "full"),
                 up = 1000,
                 down = 1000){
  
  relative = match.arg(relative)
  keytype = match.arg(keytype)
  if (keytype == 'SYMBOL' && !inherits(db_object, "OrganismDb")){
    stop("keytype of 'SYMBOL' requires input db_object to be OrganismDb object")
  }
  range_df <- OrganismDbi::select(db_object, 
                                  keys=symbol, 
                                  keytype=keytype,
                                  column = c("TXCHROM", "TXEND", "TXID", 
                                             "TXNAME", "TXSTART", "TXSTRAND"))
  
  
  if (relative == "TSS"){
    if (range_df$TXSTRAND[1] == "-"){
      start_range <- min(range_df$TXEND - down)
      end_range <- max(range_df$TXEND + up)
      strand_range <- "-"
      rel <- range_df$TXEND
    } else{
      start_range <- min(range_df$TXSTART - up)
      end_range <- max(range_df$TXSTART + down)
      strand_range <- "+"
      rel <- range_df$TXSTART
    }
    
    out <- make_view_range(chrom = range_df$TXCHROM,
                           start = start_range,
                           end = end_range,
                           strand = strand_range,
                           name = symbol,
                           relative)  
  } else if (relative == "TTS"){
    if (range_df$TXSTRAND[1] == "-"){
      start_range <- min(range_df$TXSTART - down)
      end_range <- max(range_df$TXSTART + up)
      strand_range <- "-"
      rel <- range_df$TXSTART
    } else{
      start_range <- min(range_df$TXEND - up)
      end_range <- max(range_df$TXEND + down)
      strand_range <- "+"
      rel <- range_df$TXEND
    }
    
    out <- make_view_range(chrom = range_df$TXCHROM,
                           start = start_range,
                           end = end_range,
                           strand = strand_range,
                           name = symbol,
                           relative = rel)  
    
  } else{
    if (range_df$TXSTRAND[1] == "-"){
      start_range <- min(range_df$TXSTART - down)
      end_range <- max(range_df$TXSTART + up)
      strand_range <- "-"
    } else{
      start_range <- min(range_df$TXEND - up)
      end_range <- max(range_df$TXEND + down)
      strand_range <- "+"
    }
    
    out <- make_view_range(chrom = range_df$TXCHROM,
                           start = start_range,
                           end = end_range,
                           strand = strand_range,
                           name = symbol) 
  }
  
  return(out)
}

setMethod(relative_position, signature = "ViewRange",
          function(view, position){
            position
          })

setMethod(relative_position, signature = "RelativeViewRange",
          function(view, position){
            position - mcols(view)$relative
          })







