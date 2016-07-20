#' .. content for the "Description" section of the documentation ..
#' Make data.frame that specfies where necessary files are from NGS pipeline
#' .. content for "Details" section of the documentation ..
#' @title get_chip_file_info
#'
#' @param samples_file string
#' @param pairs_file string
#'
#' @return list
#' @import gChipseq
#' @author Justin Finkle
#' @export
#' @examples
get_chip_file_info <- function(samples_file, pairs_file){
  # Load, expand, and combine sample and pair tableb
  samples <- read.delim(samples_file, stringsAsFactors = FALSE)
  samples_exp <- gChipseq::expandSampleTable(samples)
  
  pairs <- gChipseq::expandPairedTable(read.delim(pairs_file,
                                        stringsAsFactors = FALSE))
  
  full_info <- merge(pairs, samples_exp, all = TRUE)
  full_info[['Genotype']] <- samples[['Genotype']]
  
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
get_sample_names <- function(file_info){
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
add_short_names <- function(file_info, name.vector){
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
add_expression_identifier <- function(file_info, es_identifier){
file_info[['ExpressSet']] <- es_identifier
  return(file_info)
}

#' Get expression set object. Either load from rds or supply expression plot 
#' project that can be retrieved
#'
#' @param filename string: preferred method -- provide rds for ExpressionSet object. 
#' Can also provide the project name to retrieve from ExpressionPlot
#' @param stat 
#'
#' @return
#' @export
#'
#' @examples
get_expression_set <- function(filename, stat=NULL){
  # Try to load the RDS
  if(file.exists(filename)){
    eset <- readRDS(filename)
  } else if (!is.null(stat)){
    if (!requireNamespace("ExpressionPlot", quietly = TRUE)) {
      stop("ExpressionPlot is required to load a project. Please install it or load an ExpressionSet from RDS",
           call. = FALSE)
    }
    proj = ep.find.project(filename)
    eset<- ep.ExpressionSet(proj = proj, feature.type = "gene", stat = stat,
                            attach.annot = TRUE)
  } else {
    stop("No valid rds file or RNA-seq stat specified")
  }
  return(eset)
}

get_control_sample <- function(sample, file_info){
  input_name <- file_info[file_info[["Sample.Name"]] == sample, "Sample.Control"]
  return(input_name)
}

check_files <- function(df){
  sapply(grep("^File", colnames(df), value = TRUE),
         function(fileColumn) file.exists(df[, fileColumn])
  )
}

get_sample_info <- function( sample, file_info){
  sample_info <- file_info[file_info[["Sample.Name"]] == sample, ]
  return(as.list(sample_info))
}


#' Get Scaling Factor
#' Get the scaling factor for samples provided. This is assumed to be the number
#' of uniquely mapped reads
#'
#' @param file_info 
#'
#' @return
#' @export
#'
#' @examples
get_scaling_factor <- function(file_info){
  readStat1 <- ngsPipelineResult(file_info$Dir.ngs_pipeline,
                                 file='.summary_preprocess.tab$',ncores=8)
  readStat2 <- ngsPipelineResult(file_info$Dir.ngs_pipeline,ncores=8)
  readStat <- rbind(readStat1,readStat2)
  readStat.simple <- readStat[c("total_reads","analyzed"),]
  rownames(readStat.simple) <- c("Total","Uniquely mapped")
  return(readStat.simple)
}
