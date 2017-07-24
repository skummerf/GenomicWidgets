#' Load SNP data
#'
#' @description Function to load SNP information from GWAS catalog v1.0
#' @param file.name str: the file containig the snp information
#' @param sep 
#'
#' @return snp.gr GRanges: 
#' @export
#'
#' @examples
load_snp_data <- function(file.name, sep='\t'){
  snp.dat <- read.delim(file.name, sep=sep, stringsAsFactors = FALSE)
  # Remove rows that don't have a chromosome ID
  snp.dat <- snp.dat[snp.dat$CHR_ID!="",]
  snp.sqn <- unique(snp.dat$CHR_ID)
  snp.gr <- GRanges(seqnames = Rle(snp.dat$CHR_ID),
                    ranges=IRanges(start=snp.dat$CHR_POS, width=1, 
                                   names=snp.dat$SNPS))
  mcols(snp.gr) <- snp.dat
  return(snp.gr)
}