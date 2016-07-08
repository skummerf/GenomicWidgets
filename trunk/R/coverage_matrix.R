
bin_mat <- function(mat, binsize){
  tmp =  ncol(mat) / binsize
  tmp_mat = matrix(unlist(lapply(1:tmp, 
                                 function(x) c(rep(0,(x -1)*binsize),rep(1/binsize,binsize),rep(0,(tmp-x)*binsize)))),
                          nrow = ncol(mat), ncol = tmp)
  return(mat %*% tmp_mat)
}

coverage_mat_from_bigwig <- function(bigwig_file, ranges, binsize){
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
    return(tmp_mat)
  } else{
    return(bin_mat(tmp_mat, binsize))
  }
}


### This does not do read extension -- to do that you would need to use bamProfile and then smooth with
### flat window.  
coverage_mat_from_bam <- function(bam_file, ranges, binsize, ...){
  stopifnot(sum(width(ranges) == width(ranges[1])) == length(ranges)) #check for equal widths
  cvg <- bamsignals::bamCoverage(bam_file, ranges)
  tmp_mat <- t(bamsignals::alignSignals(cvg))
  if (binsize == 1){
    return(tmp_mat)
  } else{
    return(bin_mat(tmp_mat, binsize))
  }
}

coverage_mat_from_rle <- function(rle, ranges, binsize, ...){
  stopifnot(sum(width(ranges) == width(ranges[1])) == length(ranges)) #check for equal widths
  cvg <- 
    
  if (binsize == 1){
    return(tmp_mat)
  } else{
    return(bin_mat(tmp_mat, binsize))
  }
}


make_coverage_matrix <- function(inputs, 
                                 ranges, 
                                 binsize = 1, 
                                 format = c("auto","bigwig","bam", "rle","RData"),
                                 ...){
  
  format = match.arg(format)
  stopifnot(sum(width(ranges) == width(ranges[1])) == length(ranges)) #check for equal widths
  if (binsize > 1 && width(ranges[1]) %% binsize != 0){
    ranges <- resize(ranges, fix = "center", width = (width(ranges[1]) %/% binsize +1)*binsize)
  }
  
  #Determine format
  
  if (format == "bigwig"){
    # bw input
    if (length(inputs) == 1){
      out <- coverage_mat_from_bigwig(inputs, ranges, binsize)
    } else{
      out <- BiocParallel::bplapply(inputs, coverage_mat_from_bigwig, ranges, binsize)
    }
  } else if (format == "bam"){
    # bam input
    if (length(inputs) == 1){
      out <- coverage_mat_from_bam(inputs, ranges, binsize)
    } else{
      out <- BiocParallel::bplapply(inputs, coverage_mat_from_bam, ranges, binsize)
    }
  } else if (format == "RData"){
    out <- getCvgList(ranges, inputs, flank = width(ranges[1])%/%2, normalize.method = "none", bins = binsize)
  } else if (format == "rle"){
    # rle input 
    stop("not yet implemented")
  }
  return(out)
}


#Is this just the same as matrix, but with one?
get_coverage_track <- function(inputs, binsize){
  
  # bw input
  
  # bam input
  
  # rdata input
  
  
}
