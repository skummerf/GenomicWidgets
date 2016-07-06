#' .. content for the "Description" section of the documentation ..
#'
#' .. content for "Details" section of the documentation ..
#' @title getChipFileInfo
#' @param sampl_name string
#' @param samples_file string
#' @param pairs_file string
#' @return list
#' @export
#' @author Justin Finkle
#' @examples
getChipFileInfo <- function(samples_file, pairs_file){
  # Load, expand, and combine sample and pair tableb
  samples <- expandSampleTable(read.delim(samples_file, 
                                          stringsAsFactors = FALSE))
  
  pairs <- expandPairedTable(read.delim(pairs_file,
                                        stringsAsFactors = FALSE))
  
  full_info <- merge(pairs, samples, all = TRUE)
  
  return(full_info)
}


#' Title
#'
#' @param narrowpeak_file 
#' @author Justin Finkle <finklej@gene.com>
#' @return GRange
#' @export
#'
#' @examples
loadNarrowPeaks <- function(macs_peaks) {
  if (!requireNamespace('gChipseq', quietly = TRUE)){
    stop("gChipseq is required for this function to work properly. Please 
         install it", call. = FALSE)
  }
  gr_narrowpeak <- readMacsPeaks(macs_peaks)
  chr <- as.character(unique(seqnames(gr_narrowpeak)))
  
}