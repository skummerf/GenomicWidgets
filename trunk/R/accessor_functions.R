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