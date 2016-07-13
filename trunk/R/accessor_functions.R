#' .. content for the "Description" section of the documentation ..
#'
#' .. content for "Details" section of the documentation ..
#' @title getChipFileInfo
#' @param sampl_name string
#' @param samples_file string
#' @param pairs_file string
#' @return list
#' @import gChipseq
#' @author Justin Finkle
#' @export
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
#' @param file_info 
#'
#' @return
#' @export
#'
#' @examples
getSampleNames <- function(file_info){
  return(file_info[["Sample.Name"]])
}

#' Title
#'
#' @param file_info 
#' @param name.vector 
#'
#' @return
#' @export
#'
#' @examples
addShortNames <- function(file_info, name.vector){
  file_info <- data.frame(append(file_info, list("Short.Name"=name.vector), 
                                 after=match("Sample.Name", names(file_info))))
  return(file_info)
}

#' Title
#'
#' @param file_info 
#' @param es_identifier 
#'
#' @return
#' @export
#'
#' @examples
addExpressionIdentifier <- function(file_info, es_identifier){
file_info[['ExpressSet']] <- es_identifier
  return(file_info)
}

#' Title
#'
#' @param proj 
#' @import ExpressionPlot
#' @return
#' @export
#'
#' @examples
getExpressionSet <- function(filename, stat=NULL){
  # Try to load the RDS
  if(file.exists(filename)){
    eset <- readRDS(filename)
  } else if (!is.null(stat)){
    proj = ep.find.project(filename)
    eset<- ep.ExpressionSet(proj = proj, feature.type = "gene", stat = stat,
                            attach.annot = TRUE)
  } else {
    stop("No valid rds file or RNA-seq stat specified")
  }
  return(eset)
}

getControlSample <- function(sample, file_info){
  input_name <- file_info[file_info[["Sample.Name"]] == sample, "Sample.Control"]
  return(input_name)
}

checkFiles <- function(df){
  sapply(grep("^File", colnames(df), value = TRUE),
         function(fileColumn) file.exists(df[, fileColumn])
  )
}

getSampleInfo <- function( sample, file_info){
  sample_info <- file_info[file_info[["Sample.Name"]] == sample, ]
  return(as.list(sample_info))
}

