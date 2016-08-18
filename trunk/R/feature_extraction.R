#' Title
#'
#' @param gr GRanges to find peaks in
#' @param peaks_gr GRanges of Macs peaks, e.g. from readMacsPeaks. mcols must include summit and pileup
#'
#' @return
#' @export
#'
#' @examples
avg_peak_in_range <- function(gr, peaks_gr){
  peaks_summit <- GRanges(seqnames = seqnames(peaks_gr),
                          ranges = IRanges(mcols(peaks_gr)$summit, width=1),
                          score = mcols(peaks_gr)$pileup)
  overlaps <- findOverlaps(peaks_summit, gr)
  peaks_per_range <- splitAsList(mcols(peaks_summit)$score[queryHits(overlaps)],
                                 factor(subjectHits(overlaps)))
  
  # Need a better way to scale the data
  range_score <- lapply(peaks_per_range, mean)
  mcols(gr)[['avg_peak']] <- unlist(range_score)
  return(gr)
}
  


#' Title
#'
#' @param gr 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
log_scale_range <- function(gr, n = 10){
  # Split the range in two
  split_range <- tile(x = gr, n =2)[[1]]
  
  # Tile each half recursively
  # The left range is recurisively tiled on the right side, so that the granges
  left_range <- recursive_tile_range(split_range[1], side = "right", n = n)
  right_range <- recursive_tile_range(split_range[2], side = "left", n = n)
  return(c(left_range, right_range))
}

recursive_tile_range <- function(gr, 
                            side = c('right', 'left'), 
                            n = 10){
  if(n <=1){
    stop("invalid n for recursive tiling. n must be > 1")
  }
  side <- match.arg(side)
  # If the range isn't wide enough to split into n bins, you're at the bottom
  if(width(gr)<n){
    return(gr)
  } else {
    # Split into n bins
    tiled_range <- tile(x = gr, n = n )[[1]]
    
    # Call recursion on appropriate side
    if(side == 'right'){
      left_range <- tiled_range[1:(n-1)]
      right_range <- recursive_tile_range(tiled_range[n], side = side, n = n)
      
    } else if(side == 'left'){
      right_range <- tiled_range[2:n]
      left_range <- recursive_tile_range(tiled_range[1], side = side, n = n)
    }
    return(c(left_range, right_range))
  }
}