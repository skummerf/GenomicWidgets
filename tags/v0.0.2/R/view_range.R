#' Title
#'
#' @param tx_data 
#' @param db_object 
#' @param symbol 
#' @param chr 
#' @param start 
#' @param end 
#' @param strand 
#'
#' @return
#' @export
#'
#' @examples
get_view_range <- function(db_object, 
                           symbol = NULL,
                           keytype = c('SYMBOL', 'TXNAME'),
                           chr = NULL, 
                           start = NULL, 
                           end = NULL, 
                           strand = c("*", "+", "-")){
  if(is.null(symbol) && length(c(chr, start, end))<3){
    stop("Input Error. Must supply either a symbol or chr, start, and end")
  }
  strand <- match.arg(strand)
  keytype <- match.arg(keytype)
  # Make the view range from the input
  if(is.null(symbol)){
    view_range <- GRanges(Rle(c(chr)), IRanges(start = start, end = end), 
                          strand = strand)
  } else {
    view_range <- get_symbol_range(db_object = db_object, symbol = symbol, 
                                   keytype = keytype)
  }
  return(view_range)
}

#' Get Symbol Range
#' 
#' @description Make a GRange that cover the exon list
#' @param exon.df data.frame: the exons included in the annotation
#'
#' @author Justin Finkle
#' @return range.gene GRange: s
#' @export
#'
#' @examples
get_symbol_range <- function(db_object, symbol, keytype){
  db_class <- get_db_class(db_object)
  if(db_class != "OrganismDb" && keytype == "SYMBOL"){
    stop("keytype 'SYMBOL' requires db_object be of class 'OrganismDb'")
  }
  range_df <- OrganismDbi::select(db_object, keys=symbol, keytype=keytype,
                                  column = c("TXCHROM", "TXEND", "TXID", 
                                             "TXNAME", "TXSTART", "TXSTRAND"))
  # Return the range as a reduced GRanges object
  return(GenomicRanges::reduce(GRanges(range_df)))
}

#' Get DB Object
#'
#' @param db_object 
#'
#' @return
#' @export
#'
#' @examples
get_db_class <- function(db_object){
  allowed_classes <- c("TxDb", "OrganismDb")
  db_class <- class(db_object)
  if(db_class %in% allowed_classes){
    return(db_class)
  } else {
    stop(paste0("Invalid database object. Please use one of: ", 
                paste(allowed_classes, collapse = ", ")))
  }
}

#' Extend GRange
#'
#' @description extend the GRange object in either direction
#' @param gr GRange
#' @param extend scalar: 1 value if symmetric extension. Otherwise in order start, end
#'
#' @return gr GRange
#' @export
#' @author Justin Finkle
#'
#' @examples
extend_grange <- function(gr, extend=0){
  if (length(extend) == 1) {
    extend <- rep(extend,2)
  } else if ( length(extend) > 2) {
    warning("extend should be a vector of length 1 or 2. extension set to 0.")
    extend <- rep(0,2)
  }
  
  # Add the extension to the gene range for each strand
  if (as.character(strand(gr)) == '+' | as.character(strand(gr)) == '*') {
    start(gr) <- start(gr) - extend[1]
    end(gr) <- end(gr) + extend[2]
  } else if ( as.character(strand(gr)) == '-' ) {
    end(gr) <- end(gr) + extend[1]
    start(gr) <- start(gr) - extend[2]
  }
  
  return(gr)
}

#' Get Coverage in Range
#'
#' @param bwList vector: list of files to use to get coverage. Currently only accepts BigWig files
#' @param target_range GRange: specifies the range in which to get coverage
#' @param names vector: names to give each GRange coverage created
#' @param cvg_scaling vector: values by which to scale each coverage value. See get_cvg_scaling for more information.
#'
#' @return cvg_list GRangesList: coverage for each file supplied
#' @export
#' @import rtracklayer
#'
#' @examples
get_coverage_in_range <- function(bwList, 
                                  target_range, 
                                  sample_names = NULL, 
                                  cvg_scaling=NULL){
  # Import only the range that matches target_range
  cvg_list <- lapply(bwList, function(x) import.bw(x, which=target_range))
  cvg_list <- GRangesList(cvg_list)
  if(!is.null(cvg_scaling)){
    for(g in 1:length(cvg_list)){
      mcols(cvg_list[[g]])[['score']] <- mcols(cvg_list[[g]])[['score']]/cvg_scaling[[g]]
    }
  }
  
  # Name the list
  if (!missing(sample_names) && !is.null(sample_names)){
    names(cvg_list) <- sample_names
  } else {
    warning("Coverage file list supplied without names. Sequential numbering will be used")
    names(cvg_list) <- seq_along(cvg_list)
  }
  return(cvg_list)
}

get_midpoint <- function(gr){
  start(gr) + width(gr) %/% 2
}