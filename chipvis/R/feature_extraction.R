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
#' @param center 
#' @param chromosome 
#' @param n 
#' @param extension 
#' @param min_bin 
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
#'
#' @return
#' @export
#'
#' @examples
read_macs_peaks <- function(sample_list, file_info){
  peak_list <- mclapply(seq_along(sample_list), function(i, file_info){
    # The name of each treatment is being amended, so an index is passed
    # for reference
    x <- sample_list[[i]]
    new_names <- paste0(rep(names(sample_list)[[i]], length(x)), seq_along(x))
    tmp <- list()
    for(n in seq_along(x)){
      tmp[new_names[[n]]] <- readMacsPeaks(file_info$File.macs[x[[n]]])
    }
    return(tmp)
  }, file_info, mc.cores=8)
  
  # Unlist the peaks to turn it into a named list
  return(unlist(peak_list))
}

#' Title
#'
#' @param peak_list 
#' @param gr 
#' @importFrom biovizBase mold
#' @return
#' @export
#'
#' @examples
peaks_per_sample <- function(peak_list, gr){
  tmp <- gr
  for(sample in names(peak_list)){
    sample_peaks <- peak_list[[sample]]
    peaks_in_range <- avg_peak_in_range(gr = gr, peaks_gr = sample_peaks)
    mcols(tmp)[[sample]] <- mcols(peaks_in_range)[['avg_peak']]  
  }
  return(tmp)
}


# This is a terrible function. It should probably not be exported without a
# SIGNIFICANT refactor
learning <- function(s, 
                     T47D_eset,
                     truniq,
                     peaks,
                     n_recurse_tiles = 10,
                     min_bin = 100,
                     conditions = c("E2", "P4", "R5020")){
  gene <- fData(T47D_eset)$symbol[[s]]
  # print(gene)
  # Handle when genes are not in the lookup table
  if(is.na(gene) | !gene %in% truniq$symbol){
    #Need a better way to initialize zeros
    return(c(rep(0, 74)))
  }
  # Get the range
  gene_idx <- which(truniq$symbol == gene)
  gene_range <-get_view_range(chr = truniq$chr[gene_idx],
                              start = truniq$start[gene_idx],
                              end = truniq$end[gene_idx],
                              strand = as.character(truniq$strand[gene_idx]))
  pseudo_tss <- ifelse(strand(gene_range)=="-", 
                       end(gene_range), 
                       start(gene_range))
  scaled_gene_range <- log_scale_range(center = pseudo_tss, 
                                       chromosome = seqnames(gene_range),
                                       n = n_recurse_tiles,
                                       min_bin = min_bin)
  # Get the peaks
  log_peaks <- peaks_per_sample(peak_list = peaks, gr = scaled_gene_range)
  df <- biovizBase::mold(log_peaks)
  hm_vals <- apply(t(as.matrix(as.data.frame(mcols(log_peaks)))), 2, rev)+1
  
  # Get the expression
  symbol_exprs <- v$E[which(fData(T47D_eset)$symbol == gene),]
  names(symbol_exprs) <- pData(T47D_eset)$title
  symbol_exprs <- as.data.frame(symbol_exprs)
  symbol_exprs$cond <- gsub("_Rep([0-9]{1})", "", rownames(symbol_exprs))
  by_cond <- group_by(symbol_exprs, cond)
  exprs_stats <- summarize(by_cond, mean=mean(symbol_exprs), sd = sd(symbol_exprs))
  exprs_stats$cond <- conditions
  cvg_df <- as.data.frame(hm_vals)
  cvg_df$cond <- "E2"
  cvg_df$cond[grep("Pr", rownames(cvg_df))] <- "P4"
  cvg_df$cond[grep("R5020", rownames(cvg_df))] <- "R5020"
  cvg_by_cond <- group_by(cvg_df, cond)
  if(exists("filtered")){
    rm(filtered)
  }
  for(c in conditions){
    data <- rnorm(n = nrow(filter(cvg_by_cond, cond==c)),
                  mean = filter(exprs_stats, cond==c)$mean,
                  sd = filter(exprs_stats, cond==c)$sd)
    cur_sample <- filter(cvg_by_cond, cond==c)
    cur_sample$exprs <- data
    if(exists("filtered")){
      filtered <- rbind(filtered, cur_sample)
    } else{
      filtered <- cur_sample
    }
  }
  # Learn
  x <- (as.matrix(filtered[, 1:(ncol(filtered)-2)]))
  scale <- cumsum(sort(width(scaled_range)[1:37]))
  colnames(x) <- c(-rev(scale), scale)
  y <-  filtered$exprs
  result <- tryCatch({
    cvfit <- cv.glmnet(x, y)
    coefs <- coef(cvfit, s = "lambda.min")
    matrix(coefs)[2:length(coefs)]
  }, warning = function(w) {
    c(rep(0, 74))
  }, error = function(e) {
    c(rep(0, 74))
  })
  return(result)
}

