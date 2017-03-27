#' make_coverage_matrix
#' 
#' @param inputs filenames of bigwig or bam
#' @param ranges ranges for which to compute coverage within
#' @param binsize binsize to bin coverage
#' @param format format of files, default is auto 
#' @param up basepairs upstream of center to use
#' @param down basepairs downstream of center to use
#' @details up and down are 0 by default -- if not specified, actual range is used.  All ranges
#' must be of equal width.  If up and/or down are provided, then the center of the range and up 
#' basepairs upstream and down basepairs downstream are used.
#' @import GenomicRanges
#' @export
#' @author Alicia Schep
make_coverage_matrix <- function(inputs, 
                                 ranges, 
                                 binsize = 1, 
                                 format = c("auto","bigwig","bam"),
                                 up = 0,
                                 down = 0,
                                 ...){
  
  format = match.arg(format)
  if (up > 0 || down > 0){
    ranges <- resize(ranges, fix = "center", width = 1)
    ranges <- promoters(ranges, upstream = up, downstream = down)
  } 
  
  stopifnot(sum(width(ranges) == width(ranges[1])) == length(ranges)) #check for equal widths
  if (binsize > 1 && width(ranges[1]) %% binsize != 0){
    ranges <- resize(ranges, fix = "center", width = (width(ranges[1]) %/% binsize +1)*binsize)
  }
  
  if (up > 0 || down > 0){
    coln <- seq(-up,down - binsize,binsize)
  } else{
    coln <- seq(1,width(ranges[1]),binsize)
  }
  
  #Determine format
  if (format == "auto"){
    if (length(inputs) == 1){
      tmp = inputs
    } else{
      tmp = inputs[[1]]
    }
    if (is.character(tmp)){
      if (substr(tmp, nchar(tmp) - 6 , nchar(tmp)) == ".bigwig" || 
          substr(tmp, nchar(tmp) - 2 , nchar(tmp)) == ".bw"){
        format = "bigwig"
      } else if (substr(tmp, nchar(tmp) - 3 , nchar(tmp)) == ".bam"){
        format = "bam"
      } else{
        stop("Cannot determine format of inputs.")
      }
    }
  }
  
  if (format == "bigwig"){
    # bw input
    if (length(inputs) == 1){
      out <- coverage_mat_from_bigwig(inputs, ranges, binsize, coln)
    } else{
      out <- lapply(inputs, coverage_mat_from_bigwig, ranges, binsize, coln)
    }
  } else if (format == "bam"){
    # bam input
    #inputs = convert_to_full_path(inputs) #bamsignals seems to fail with symlinks
    #
    if (length(inputs) == 1){
      out <- coverage_mat_from_bam(inputs, ranges, binsize, coln)
    } else{
      out <- lapply(inputs, coverage_mat_from_bam, ranges, binsize, coln)
    }
  } else {
    stop("Format not recognized")
  }
  
  return(out)
}

#' normalize_coverage_matrix
#' 
#' @param mats
#' @param method normalization method option, see Details
#' @param pct Percentile, only used if PercentileMax is method
#' @param scalar vector of scalars used for normalizing each mat, only 
#' used if scalar is method
#' @details Normalization choices are "localRms", "localMean", 
#' "localNonZeroMean", "PercentileMax", "scalar", and "none".  
#' @export
#' @author Alicia Schep
normalize_coverage_matrix <- function(mats, 
                                      method = c("localRms", "localMean", 
                                                      "localNonZeroMean", "PercentileMax", "scalar", "none"), 
                                      pct = 0.95, scalar = NULL){
  method = match.arg(method)
  if (!is.list(mats)){
    out <- normalize_coverage_matrix_single(mats, method, pct, scalar)
  } else{
    if (method == "scalar"){
      stopifnot(!is.null(scalar))
        out <- lapply(1:length(mats), function(x) {
          normalize_coverage_matrix_single(mats[[x]],method = "scalar",scalar = scalar[x])
        })
        names(out) <- names(mats)
    } else{
      out <- lapply(mats, normalize_coverage_matrix_single, method = method, pct = pct)
    }
  }
  return(out)
}


## Modified from gChipseq package
## not exported
normalize_coverage_matrix_single <- function(mat, method = c("localRms", "localMean", 
                                     "localNonZeroMean", "PercentileMax", "scalar", "none"), 
          pct = 0.95, scalar = NULL) 
{
  method <- match.arg(method)
  if (method == "localRms") {
    scaler <- sqrt(rowSums(mat^2))
    scaler[scaler == 0] <- 1
    return(mat/scaler)
  } else if (method == "localMean") {
    scaler <- rowMeans(mat)
    scaler[scaler == 0] <- 1
    return(mat/scaler)
  } else if (method == "localNonZeroMean") {
    scaler <- apply(mat, 1, function(x) {
      mean(x[x != 0])
    })
    scaler[is.nan(scaler)] <- 1
    return(mat/scaler)
  } else if (method == "PercentileMax") {
    scaler <- quantile(mat, pct)
    if (scaler == 0) {
      scaler <- 1
    }
    res <- mat/scaler
    res[res[] > 1] <- 1
    return(res)
  } else if (method == "scalar") {
    return(mat/scalar)
  } else if (method == "none") {
    return(mat)
  }
}



#### Helper Functions (not exported) ---------------------------------------------------------

#Function to reduce size of matrix by binning columns
bin_mat <- function(mat, binsize){
  out <- sapply(seq(1, ncol(mat),binsize), function(x) rowMeans(mat[,x:(x+binsize - 1),drop =FALSE]))
  if (dim(mat)[1] == 1) out <- matrix(out, nrow = 1)
  return(out)
}


coverage_mat_from_bigwig <- function(bigwig_file, ranges, binsize, coln){
  stopifnot(sum(width(ranges) == width(ranges[1])) == length(ranges)) #check for equal widths
  tmp = rtracklayer::import.bw(bigwig_file, selection = ranges, as = "NumericList")
  tmp_list = lapply(1:length(ranges), function(x){
    if (as.character(strand(ranges[x])) == "-"){
      return(rev(tmp@listData[[x]]))
    } else{
      return(tmp@listData[[x]])
    }
  }) 
  tmp_mat = t(simplify2array(tmp_list))
  if (binsize == 1){
    colnames(tmp_mat) = coln
    return(tmp_mat)
  } else{
    binned <- bin_mat(tmp_mat, binsize)
    colnames(binned) <- coln
    return(binned)
  }
}


### This does not do read extension -- to do that you would need to use bamProfile and then smooth with
### flat window.  
coverage_mat_from_bam <- function(bam_file, ranges, binsize, coln){
  stopifnot(sum(width(ranges) == width(ranges[1])) == length(ranges)) #check for equal widths
  tmp_mat <- bamsignals::bamCoverage(bam_file, ranges, verbose = FALSE) %>% bamsignals::alignSignals() %>% t()
  #tmp_mat <- t(bamsignals::alignSignals(cvg))
  if (binsize == 1){
    colnames(tmp_mat) = coln
    return(tmp_mat)
  } else{
    binned <- bin_mat(tmp_mat, binsize)
    colnames(binned) <- coln
    return(binned)
  }
}


#Function from Johannes Helmuth, https://github.com/lamortenera/bamsignals/issues/16
bamCoverageFrags <- function(bampath, gr, fragLength=200, ss = TRUE, ...) {
  prf <- bamsignals::bamProfile(bampath, gr, binsize=1, ss = TRUE, ...)
  cov <- lapply(as.list(prf), function(p) {
    l <- dim(p)[2]
    x1 <- cumsum(as.numeric(p[1,]))
    y1 <- c(rep(0, fragLength), x1[1:(l-200)])
    x2 <- cumsum(as.numeric(rev(p[2,])))
    y2 <- c(rep(0, fragLength), x2[1:(l-200)])
    matrix(c(x1-y1, rev(x2-y2)), ncol = l, byrow = T,
           dimnames = list(c("sense", "antisense"), NULL))
  })
  if (ss == TRUE) {
    prf@signals <- cov
  } else {
    prf@signals <- lapply(cov, colSums)
  }
  invisible(prf)
}


convert_to_full_path <- function(x){
  for (i in 1:length(x)){
    x[[i]] = system(paste0("readlink -f ", x[[i]]), intern = TRUE)
  }
  x
}

bin_track_mat <- function(track_mat, target_range, tiled_range, scaling_factors){
  w <- width(tiled_range)
  wcs <- c(0,cumsum(w))
  out <- sapply(seq_along(tiled_range), function(x){
    s <- wcs[x] + 1
    e <- wcs[x + 1]
    colMeans(track_mat[s:e,, drop = FALSE], na.rm = TRUE)
  } )
  out <- t(out / scaling_factors)
  return(out)
}


coverage_track_from_bigwig <- function(bigwig_file, target_range){
  rtracklayer::import.bw(bigwig_file, selection = target_range, as = "NumericList")@listData[[1]]
}

coverage_track_from_bam <- function(bam_file, target_range){
  as.list(bamsignals::bamCoverage(bam_file, target_range, verbose = FALSE))[[1]]
}



#' make_coverage_tracks
#' @param inputs filenames of bigwig, bam, or RData file
#' @param ranges ranges for which to compute coverage within
#' @param binsize binsize to bin coverage
#' @param format format of files, default is auto 
#' @param up basepairs upstream of center to use
#' @param down basepairs downstream of center to use
#' @details up and down are 0 by default -- if not specified, actual range is used.  All ranges
#' must be of equal width.  If up and/or down are provided, then the center of the range and up 
#' basepairs upstream and down basepairs downstream are used.
#' @import GenomicRanges
#' @export
#' @author Alicia Schep
make_coverage_tracks <- function(inputs, 
                                 target_range, 
                                 sample_names = names(inputs),
                               binsize = 100, 
                               method = c("a","j"),
                               format = c("auto","bigwig","bam"), 
                               bin = TRUE,
                               scaling_factors = rep(1, length(inputs)),
                               up = 0,
                               down = 0,
                               stranded = FALSE){
  format = match.arg(format)
  method = match.arg(method)
  stopifnot(length(target_range) == 1)
  negative_strand <- FALSE
  
  if (up > 0 || down > 0){
      target_range <- resize(target_range, fix = "center", width = 1)
      target_range <- promoters(target_range, upstream = up, downstream = down)
    } 
  #if (binsize > 1 && width(target_range[1]) %% binsize != 0){
  #    target_range <- resize(target_range, fix = "center", width = (width(target_range[1]) %/% binsize +1)*binsize)
  #}
  if (stranded && strand(target_range) == "-"){
    negative_strand <- TRUE
  }
  if (format == "auto"){
    if (length(inputs) == 1){
      tmp = inputs
    } else{
      tmp = inputs[[1]]
    }
    if (is.character(tmp)){
      if (substr(tmp, nchar(tmp) - 6 , nchar(tmp)) == ".bigwig" || 
          substr(tmp, nchar(tmp) - 2 , nchar(tmp)) == ".bw"){
        format = "bigwig"
      } else if (substr(tmp, nchar(tmp) - 3 , nchar(tmp)) == ".bam"){
        format = "bam"
      } else{
        stop("Cannot determine format of inputs.")
      }
    }
  }
  strand(target_range) <- "*"
  if (method == "a"){
    if (format == "bigwig"){
      tracks <- sapply(inputs, coverage_track_from_bigwig, target_range)
    } else if (format == "bam"){
      tracks <- sapply(inputs, coverage_track_from_bam, target_range)
    } else{
      stop("Incorrect format!")
    }
    if (is.null(dim(tracks))) tracks <- matrix(tracks, nrow = 1)
    colnames(tracks) <- sample_names
    out <- GenomicRanges::tile(target_range, width = binsize)[[1]]
    tracks <- bin_track_mat(tracks, target_range, out, scaling_factors)
    mcols(out) <- tracks
    if (stranded && negative_strand) out <- rev(out)
  } else{
    stopifnot(format == "bigwig")
    out <- get_coverage_in_range(inputs, 
                                 target_range = target_range, 
                                 cvg_scaling = scaling_factors,
                                 sample_names = sample_names)
    out <- bin_coverage_in_range(out, target_range, binwidth = binsize)
    if (stranded && negative_strand) out <- rev(out)
  }
  return(out)
}

#' subset_coverage
#' 
#' @export
#' @author Alicia Schep
subset_coverage <- function(coverage, idx){
  mcols(coverage) <- mcols(coverage)[,idx]
  coverage
}
