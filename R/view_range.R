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
                 keytype = c("SYMBOL", "TXNAME"),
                 relative = c("TSS", "TTS", "full"),
                 up = 1000,
                 down = 1000){
  
  relative <- match.arg(relative)
  keytype <- match.arg(keytype)
  if (keytype == "SYMBOL" && !inherits(db_object, "OrganismDb")){
    stop("keytype of 'SYMBOL' requires input db_object to be OrganismDb object")
  }
  range_df <- OrganismDbi::select(db_object, 
                                  keys=symbol, 
                                  keytype=keytype,
                                  column = c("TXCHROM", "TXEND", "TXID", 
                                             "TXNAME", "TXSTART", "TXSTRAND"))
  
  if (length(unique(range_df$TXCHROM)) > 1){
    stop("symbol mapped to multiple chromosomes")
  }
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
    
  }
  out <- GRanges(range_df$TXCHROM[1], 
                 ranges = IRanges(
                   start = start_range,
                   end = end_range),
                 strand = strand_range,
                 name = symbol) 
  
  return(out)
}

setMethod(relative_position, signature = "ViewRange",
          function(view, position){
            stopifnot(length(view) == 1)
            if (view@relative){
              if (as.character(strand(view@range)) == "-"){
                view@reference - position
              } else{
                position - view@reference 
              }
            } else{
              position
            }
          })







