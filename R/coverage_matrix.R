#' make_coverage_matrix
#' 
#' Makes coverage matrix based on bam or bigwig inputs
#' 
#' @param inputs filenames of bigwig or bam
#' @param ranges ranges for which to compute coverage within
#' @param input_names names to associate with input bigwig or bam files, used 
#' for naming assays in resulting SummarizedExperiment object
#' @param binsize binsize to bin coverage
#' @param format format of files, default is auto 
#' @param up basepairs upstream of center to use
#' @param down basepairs downstream of center to use
#' @details up and down are 0 by default -- if not specified, actual range is
#'  used.  All ranges must be of equal width.  If up and/or down are provided,
#'  then the center of the range and up basepairs upstream and down basepairs
#'  downstream are used.
#' @return RangedSummarizedExperiment 
#' @import GenomicRanges
#' @export
#' @author Alicia Schep
#' @examples
#' 
#' library(GenomicRanges)
#' # First read in some sample data
#' genomation_dir <- system.file("extdata", package = "genomationData")
#'
#' samp.file <- file.path(genomation_dir,'SamplesInfo.txt')
#' samp.info <- read.table(samp.file, header=TRUE, sep='\t', 
#'                        stringsAsFactors = FALSE)
#' samp.info$fileName <- file.path(genomation_dir, samp.info$fileName)
#'
#' ctcf.peaks = genomation::readBroadPeak(system.file("extdata",
#'                                                   "wgEncodeBroadHistoneH1hescCtcfStdPk.broadPeak.gz",
#'                                                   package = "genomationData"))
#' ctcf.peaks = ctcf.peaks[seqnames(ctcf.peaks) == "chr21"]
#' ctcf.peaks = ctcf.peaks[order(-ctcf.peaks$signalValue)]
#' ctcf.peaks = resize(ctcf.peaks, width = 1000, fix = "center")
#'
#' # Make the coverage matrices
#' mats <- make_coverage_matrix(samp.info$fileName[1:3], ctcf.peaks, 
#'                      up = 500, down = 500, binsize = 25)
#'
#' # Benchmarking speed of make_coverage_matrix compared to ScoreMatrixList
#' function from genomation
#' \dontrun{
#' bm <- microbenchmark::microbenchmark(ctcf_mats = 
#'                              make_coverage_matrix(samp.info$fileName[1:3], 
#'                                                     ctcf.peaks, 
#'                                                    up = 500, down = 500, 
#'                                                     binsize = 25), 
#'                              geno = ScoreMatrixList(samp.info$fileName[1:3], 
#'                                               ctcf.peaks, bin.num = 1000/25),
#'                              times = 5)
#'
#' bm
#'
#' plot(bm)}
#'
make_coverage_matrix <- function(inputs, 
                                 ranges, 
                                 input_names = names(inputs),
                                 binsize = 1, 
                                 format = c("auto","bigwig","bam"),
                                 up = 0,
                                 down = 0){
  
  format <- match.arg(format)
  
  if (up > 0 || down > 0){
    ranges <- resize(ranges, fix = "center", width = 1)
    ranges <- promoters(ranges, upstream = up, downstream = down)
  } 
  #check for equal widths
  stopifnot(sum(width(ranges) == width(ranges[1])) == length(ranges)) 
  if (binsize > 1 && width(ranges[1]) %% binsize != 0){
    ranges <- resize(ranges, fix = "center", 
                     width = (width(ranges[1]) %/% binsize +1)*binsize)
  }
  
  rn <- paste(as.character(seqnames(ranges)), 
              paste(end(ranges),start(ranges), sep="-"), 
              sep=":")
  names(ranges) <- rn
  
  if (up > 0 || down > 0){
    coln <- seq(-up,down - binsize,binsize)
  } else{
    coln <- seq(1,width(ranges[1]),binsize)
  }
  
  names(inputs) <- input_names
  #Determine format
  if (format == "auto"){
    if (length(inputs) == 1){
      tmp <- inputs
    } else{
      tmp <- inputs[[1]]
    }
    if (is.character(tmp)){
      if (substr(tmp, nchar(tmp) - 6 , nchar(tmp)) == ".bigwig" || 
          substr(tmp, nchar(tmp) - 2 , nchar(tmp)) == ".bw"){
        format <- "bigwig"
      } else if (substr(tmp, nchar(tmp) - 3 , nchar(tmp)) == ".bam"){
        format <- "bam"
      } else{
        stop("Cannot determine format of inputs.")
      }
    }
  }
  
  if (format == "bigwig"){
    # bw input
    if (length(inputs) == 1){
      out <- coverage_mat_from_bigwig(inputs, ranges, binsize, coln)
      rownames(out) <- rn
      out <- list(out)
      names(out) <- input_names
    } else{
      out <- lapply(inputs, coverage_mat_from_bigwig, ranges, binsize, coln)
      out <- lapply(out, function(x) {rownames(x) <- rn; x})
    }
  } else if (format == "bam"){

    if (length(inputs) == 1){
      out <- coverage_mat_from_bam(inputs, ranges, binsize, coln)
      rownames(out) <- rn
      out <- list(out)
      names(out) <-  input_names
    } else{
      out <- lapply(inputs, coverage_mat_from_bam, ranges, binsize, coln)
      out <- lapply(out, function(x) {rownames(x) <- rn; x})
    }
  } else {
    stop("Format not recognized")
  }
  return(SummarizedExperiment::SummarizedExperiment(out, rowRanges = ranges))
}

#' normalize_coverage_matrix
#' 
#' Normalizes coverage matrices using one of several methods.
#' @param mats matrix, list of matrix, or SummarizedExperiment
#' @param method normalization method option, see Details
#' @param pct Percentile, only used if PercentileMax is method
#' @param scalar vector of scalars used for normalizing each mat, only 
#' used if scalar is method
#' @param ... additional arguments to normalize_coverage_matrix
#' @details Normalization choices are "localRms", "localMean", 
#' "localNonZeroMean", "PercentileMax", "scalar", and "none".  localRMS will 
#' divide each row by the root mean squared values of that row.  localMean will
#' divide each row by the mean of that row.  localNonZeroMean will divide each 
#' row by nonzero values in that row.  PercentileMax will divide values based on 
#' percentile (given by pct argument) of the entire matrix.  scalar will divide
#' entire matrix by a scalar, given by scalar argument.  This scalar could for 
#' example be a measure of the sequencing depth.  
#' @export
#' @rdname normalize_coverage_matrix
#' @name normalize_coverage_matrix
#' @aliases normalize_coverage_matrix,list-method
#' normalize_coverage_matrix,matrix-method
#' normalize_coverage_matrix,SummarizedExperiment-method
#' @author Alicia Schep
setMethod(normalize_coverage_matrix, "list",
          function(mats, 
                   method = c("localRms", 
                              "localMean", 
                              "localNonZeroMean", 
                              "PercentileMax", 
                              "scalar", "none"), 
                   pct = 0.95, 
                   scalar = NULL){
            method <- match.arg(method)
            if (method == "scalar"){
              stopifnot(!is.null(scalar))
              out <- lapply(seq_len(length(mats)), function(x) {
                normalize_coverage_matrix_single(mats[[x]],
                                                 method = "scalar",
                                                 scalar = scalar[x])
              })
              names(out) <- names(mats)
            } else{
              out <- lapply(mats, normalize_coverage_matrix_single, 
                            method = method, pct = pct)
            }
            
            return(out)
          })
#' @rdname normalize_coverage_matrix
#' @export
setMethod(normalize_coverage_matrix, "matrix",
          function(mats, 
                   method = c("localRms", 
                              "localMean", 
                              "localNonZeroMean", 
                              "PercentileMax", 
                              "scalar", "none"), 
                   pct = 0.95, 
                   scalar = NULL){
            method <- match.arg(method)
            out <- normalize_coverage_matrix_single(mats, method, pct, scalar)
            return(out)
          })

#' @export
#' @rdname normalize_coverage_matrix
setMethod(normalize_coverage_matrix, "SummarizedExperiment",
          function(mats, 
                   ...){
            out <- mats
            assays(out) <- 
              normalize_coverage_matrix(as.list(assays(mats)), ...)
            return(out)
          })

## Modified from gChipseq package
## not exported
normalize_coverage_matrix_single <- function(mat, 
                                             method = c("localRms", 
                                                        "localMean", 
                                                        "localNonZeroMean", 
                                                        "PercentileMax", 
                                                        "scalar", 
                                                        "none"), 
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
  out <- vapply(seq(1, ncol(mat),binsize), 
                function(x) rowMeans(mat[,x:(x+binsize - 1),drop =FALSE]),
                rep(0,nrow(mat)))
  if (dim(mat)[1] == 1) out <- matrix(out, nrow = 1)
  return(out)
}


coverage_mat_from_bigwig <- function(bigwig_file, ranges, binsize, coln){
  #check for equal widths
  stopifnot(sum(width(ranges) == width(ranges[1])) == length(ranges)) 
  tmp <- rtracklayer::import.bw(bigwig_file, selection = ranges, 
                               as = "NumericList")
  tmp_list <- lapply(seq_along(ranges), function(x){
    if (as.character(strand(ranges[x])) == "-"){
      return(rev(tmp@listData[[x]]))
    } else{
      return(tmp@listData[[x]])
    }
  }) 
  tmp_mat <- t(simplify2array(tmp_list))
  if (binsize == 1){
    colnames(tmp_mat) <- coln
    return(tmp_mat)
  } else{
    binned <- bin_mat(tmp_mat, binsize)
    colnames(binned) <- coln
    return(binned)
  }
}


# This does not do read extension -- to do that you would need to use 
# bamProfile and then smooth with flat window.  
#' @importFrom bamsignals bamCoverage alignSignals
coverage_mat_from_bam <- function(bam_file, ranges, binsize, coln){
  #check for equal widths
  stopifnot(sum(width(ranges) == width(ranges[1])) == length(ranges)) 
  tmp_mat <- t(alignSignals(bamCoverage(bam_file, ranges, verbose = FALSE)))
  #tmp_mat <- t(bamsignals::alignSignals(cvg))
  if (binsize == 1){
    colnames(tmp_mat) <- coln
    return(tmp_mat)
  } else{
    binned <- bin_mat(tmp_mat, binsize)
    colnames(binned) <- coln
    return(binned)
  }
}


# Function from Johannes Helmuth, 
# https://github.com/lamortenera/bamsignals/issues/16
bamCoverageFrags <- function(bampath, gr, fragLength=200, ss = TRUE, ...) {
  prf <- bamsignals::bamProfile(bampath, gr, binsize=1, ss = TRUE, ...)
  cov <- lapply(as.list(prf), function(p) {
    l <- dim(p)[2]
    x1 <- cumsum(as.numeric(p[1,]))
    y1 <- c(rep(0, fragLength), x1[seq_len(l-200)])
    x2 <- cumsum(as.numeric(rev(p[2,])))
    y2 <- c(rep(0, fragLength), x2[seq_len(l-200)])
    matrix(c(x1-y1, rev(x2-y2)), ncol = l, byrow = TRUE,
           dimnames = list(c("sense", "antisense"), NULL))
  })
  if (ss == TRUE) {
    prf@signals <- cov
  } else {
    prf@signals <- lapply(cov, colSums)
  }
  invisible(prf)
}





