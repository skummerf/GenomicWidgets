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
  
  # Average peaks in each bin
  range_score <- lapply(peaks_per_range, mean)
  mcols(gr)[['avg_peak']] <- 0
  mcols(gr)[['avg_peak']][as.numeric(names(range_score))] <- unlist(range_score)
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
log_scale_range <- function(center,
                            chromosome,
                            n = 10, 
                            extension = 1000000, 
                            min_bin = 100){
  
  new_gr <- GRanges(seqnames = chromosome,
                    ranges = IRanges(center-extension, width=2*extension))
  split_range <- tile(x = new_gr, n = 2)[[1]]
  
  # Tile each half recursively
  # The left range is recurisively tiled on the right side, so that the granges
  left_range <- recursive_tile_range(split_range[1], 
                                     side = "right", 
                                     n = n, 
                                     min_bin = min_bin)
  right_range <- recursive_tile_range(split_range[2], 
                                      side = "left", 
                                      n = n,
                                      min_bin = min_bin)
  return(c(left_range, right_range))
}

recursive_tile_range <- function(gr, 
                                 side = c('right', 'left'), 
                                 n = 10,
                                 min_bin = 100){
  if(n <=1){
    stop("invalid n for recursive tiling. n must be > 1")
  }
  side <- match.arg(side)
  # If the range isn't wide enough to split into n bins, you're at the bottom
  if(width(gr)<n | width(gr)< min_bin*n){
    return(gr)
  } else {
    # Split into n bins
    tiled_range <- tile(x = gr, n = n )[[1]]
    
    # Call recursion on appropriate side
    if(side == 'right'){
      left_range <- tiled_range[1:(n-1)]
      right_range <- recursive_tile_range(tiled_range[n], 
                                          side = side, 
                                          n = n,
                                          min_bin = min_bin)
      
    } else if(side == 'left'){
      right_range <- tiled_range[2:n]
      left_range <- recursive_tile_range(tiled_range[1], 
                                         side = side, 
                                         n = n,
                                         min_bin = min_bin)
    }
    return(c(left_range, right_range))
  }
}

#' Title
#'
#' @param symbol_range 
#' @param scaled_range 
#' @param hm_vals 
#' @param base 
#'
#' @return
#' @export
#'
#' @examples
make_lines <- function(symbol_range, scaled_range, hm_vals, base){
  n_bins <- ncol(hm_vals)
  center <- (start(symbol_range) + end(symbol_range))/2
  left_width <- floor(log(cumsum(sort(width(scaled_range[1:floor(n_bins/2)]))), base))
  l_log_scale <- unique(left_width)[unique(left_width) >2 ]
  left_lines <- sapply(l_log_scale, function(x, widths){
    first(which(widths==x))
  }, widths=left_width)
  right_width <- floor(log(cumsum(sort(width(scaled_range[(floor(n_bins/2)+1):n_bins]))), base))
  r_log_scale <- unique(left_width)[unique(left_width) >2 ]
  right_lines <- sapply(r_log_scale, function(x, widths){
    first(which(widths==x))
  }, widths=right_width)
  center_line_x <- ncol(hm_vals)/2-0.5
  binsize <- floor(log(width(scaled_range), base))
  log_lines <- sort(c(-left_lines, 0, right_lines)) + center_line_x
  names(log_lines) <-  sapply(sort(c(-l_log_scale, 0, r_log_scale)), function(x) {
    if(x > 0){
      sign_val <- "+"
    } else if( x< 0){
      sign_val <- "-"
    } else {
      sign_val <- ""
    }
    return(paste0(sign_val, base, "<sup>", abs(x), "</sup>"))
  })
  line_list <- list()
  for(ll in seq_along(log_lines)){
    x <- log_lines[[ll]]
    if(x == center_line_x){
      color <- 'red'
    } else {
      color <- 'black'
    }
    line_list[[ll]] <- list(type='line', line=list(width=1, color=color), x0=x, x1=x, y0=0, y1=1, yref='paper')
  }
  scale_list <- list(lines = line_list, vals = log_lines)
  return(scale_list)
}

## This probably shouldn't be a function. It's too specific
#' Title
#'
#' @param sample_list 
#' @param file_info 
#' @param gr 
#'
#' @return
#' @export
#'
#' @examples
peaks_per_sample <- function(sample_list, file_info, gr){
  mclapply(seq_along(sample_list), function(i, file_info, gr){
    x <- sample_list[[i]]
    new_names <- paste0(rep(names(sample_list)[[i]], length(x)), seq_along(x))
    tmp <- gr
    for(n in seq_along(x)){
      peaks <- readMacsPeaks(file_info$File.macs[x[[n]]])
      peaks_in_range <- avg_peak_in_range(gr = gr, peaks_gr = peaks)
      mcols(tmp)[[new_names[n]]] <- mcols(peaks_in_range)[['avg_peak']]
      }
    return(mcols(tmp))
    }, file_info, gr, mc.cores=8)
}