setAs("ViewRange","GRanges",
      function(from){
        out <- GRanges(seqnames = seqnames(from),
                strand = strand(from),
                ranges = IRanges(start = start(from), end = end(from)),
                seqinfo = seqinfo(from))
        mcols(out) <- mcols(from)
        return(out)
      })

setAs("GRanges","ViewRange",
      function(from){
        if (!("name" %in% colnames(mcols(from)))) 
          mcols(from)$name <- as.character(seq_along(from))
        new("ViewRange", from)
      })

setAs("GRanges","RelativeViewRange",
      function(from){
        if (!("name" %in% colnames(mcols(from)))) 
          mcols(from)$name <- ""
        stopifnot("relative" %in% colnames(mcols(from)))
        stopifnot("reference" %in% colnames(mcols(from)))
        new("RelativeViewRange", from)
      })

make_view_range <- function(chrom,
                       start, 
                       end, 
                       strand = c("*", "+", "-"),
                       relative = NULL,
                       name = character(0)){
  
  strand <- match.arg(strand)
  gr <- GRanges(Rle(c(chrom)), IRanges(start = start, end = end), 
                strand = strand)
  if (is.null(relative)){
    out <- new("ViewRange", 
               gr, 
               name = ifelse(length(name) == 0, as.character(gr), name))
  } else{
    out <- new("RelativeViewRange", 
               gr, 
               relative = relative,
               name = ifelse(length(name) == 0, as.character(gr), name))
  }
  out
}

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
    if (range_df$TXSTRAND == "-"){
      start_range <- range_df$TXEND - down
      end_range <- range_df$TXEND + up
      strand_range <- "-"
      rel <- range_df$TXEND
    } else{
      start_range <- range_df$TXSTART - up
      end_range <- range_df$TXSTART + down
      strand_range <- range_df$TXSTRAND
      rel <- range_df$TXSTART
    }
    
    out <- make_view_range(chrom = range_df$TXCHROM,
                           start = start_range,
                           end = end_range,
                           strand = strand_range,
                           name = symbol,
                           relative = rel,
                           reference = "TSS")  
  } else if (relative == "TTS"){
    if (range_df$TXSTRAND == "-"){
      start_range <- range_df$TXSTART - down
      end_range <- range_df$TXSTART + up
      strand_range <- "-"
      rel <- range_df$TXSTART
    } else{
      start_range <- range_df$TXEND - up
      end_range <- range_df$TXEND + down
      strand_range <- range_df$TXSTRAND
      rel <- range_df$TXEND
    }
    
    out <- make_view_range(chrom = range_df$TXCHROM,
                           start = start_range,
                           end = end_range,
                           strand = strand_range,
                           name = symbol,
                           relative = rel,
                           reference = "TTS")  
    
  } else{
    if (range_df$TXSTRAND == "-"){
      start_range <- range_df$TXSTART - down
      end_range <- range_df$TXSTART + up
      strand_range <- "-"
    } else{
      start_range <- range_df$TXEND - up
      end_range <- range_df$TXEND + down
      strand_range <- range_df$TXSTRAND
    }
    
    out <- make_view_range(chrom = range_df$TXCHROM,
                           start = start_range,
                           end = end_range,
                           strand = range_df$TXSTRAND,
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







