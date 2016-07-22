#' get_gene_range
#' 
#' @description Make a GRange that cover the exon list
#' @param exon.df data.frame: the exons included in the annotation
#'
#' @author Justin Finkle
#' @return range.gene GRange: s
#' @export
#'
#' @examples
get_gene_range <- function(exon.df){
  # Make the range
  range.gene <- GRanges(seqnames=exon.df[1,"chr"],
                        ranges=IRanges(start=min(exon.df$start),
                                       end=max(exon.df$end)),
                        strand=exon.df[1,"strand"])
  return(range.gene)
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
#' @param target.range GRange: specifies the range in which to get coverage
#' @param names vector: names to give each GRange coverage created
#' @param scaling.factor vector: values by which to scale each coverage value. See get_scaling_factor for more information.
#'
#' @return cvg.L GRangesList: coverage for each file supplied
#' @export
#'
#' @examples
get_coverage_in_range <- function(bwList, target.range, names, scaling.factor=NULL){
  # Import only the range that matches target.range
  cvg.L <- lapply(bwList, function(x) import.bw(x, which=target.range))
  cvg.L <- GRangesList(cvg.L)
  if(!is.null(scaling.factor)){
    for(g in 1:length(cvg.L)){
      mcols(cvg.L[[g]])[['score']] <- mcols(cvg.L[[g]])[['score']]/scaling.factor[[g]]
    }
  }
  
  # Name the list
  if (!missing(names)){
    names(cvg.L) <- names
  }
  return(cvg.L)
}

##' get exon information from IGIS
##'
##' @param src source of id, gene|symbol
##' @param org organism, human or mouse
##' @param target character or GRange: gene symbol or GRange object for which to get exon information 
##' @param genome 
##' @param ... additional argument passed to igisGetmart()
##'
##' @return a data frame of exon information
##' @import biomaRt gChipseq
##' @export
##' @author Jinfeng Liu
get_exons <- function(target, src=c('symbol','gene'), genome,
                      org=c('human','mouse'),...) {
  
  
  org <- match.arg(org)
  dataset <- ifelse(org == 'human', "hsapiens_gene_ensembl","mmusculus_gene_ensembl")
  
  attr <- c("chromosome_name","exon_chrom_start","exon_chrom_end","strand",
            "transcript_id","gene_id","external_gene_id","exon_id","rank",
            "cds_start","cds_end")
  columns <- c("chr","start","end","strand",
               "transcript","gene","symbol","exon","rank",
               "cds_start","cds_end")
  
  mart <- gChipseq::igisGetmart(...)
  mart <- biomaRt::useDataset(dataset, mart)
  
  if (is(target, "GRanges")){
    region <- list(as.integer(as.character((seqnames(target)))), 
                   min(start(target)), max(end(target)))
    exon.table <- getBM(attributes = attr,
                        filters = c("chromosome_name","start","end"),
                        values = region, mart = mart)
  } else if (is(target, "character")){
    src <- match.arg(src)
    src_filter <- ifelse(src == 'gene', 'gene_id', 'external_gene_id')
    exon.table <- getBM(attributes = attr,
                        filters = c("source",src_filter),
                        values = list("entrez",target),mart = mart)
  } else {
    warning("Must use valid entrez gene id or GRanges object")
    return(NULL)
  }
  
  colnames(exon.table) <- columns
  
  exon.table$strand <- ifelse(exon.table$strand == 1, '+',
                              ifelse(exon.table$strand == -1, '-',NA))
  
  return(exon.table)
}

#' Title
#'
#' @param symbol 
#' @param org 
#' @param genome 
#' @param chr 
#' @param start 
#' @param end 
#' @param strand 
#'
#' @return
#' @export
#'
#' @examples
get_view_range <- function(symbol, org, genome, chr=NULL, start=NULL, end=NULL, 
                           strand=NULL){
  if(!missing(symbol)){
    exon.df <- chipVis::get_exons(target=symbol, org=org, genome=genome)
    view.range <- get_gene_range(exon.df)
  } else {
    view.range <- GRanges(Rle(c(chr)), IRanges(start=start, end=end), 
                          strand = strand)
  }
  return(view.range)
}