
#' get_tss_from_txdb
#' 
#' @param txdb txdb object
#' @param by_gene return a GRangesList with tss grouped by gene
#' @param db basename of meta-data package to use for \code{\link[annotate]{getSYMBOL}}
#' @importFrom S4Vectors split
#' @return GRanges or GRangesList
#' @export
get_tss_from_txdb <- function( txdb, by_gene = TRUE, db = 'org.Hs.eg'){
  tr <- transcriptsBy(txdb, by = "gene") %>% unlist()
  tr$entrez <- stringr::str_extract(names(tr),"(?<=:)(.*)")
  tr$symbol <- annotate::getSYMBOL(tr$entrez, data=db)
  tr <- resize(tr, fix = "start", width = 1)
  names(tr) = NULL
  if (!isTRUE(by_gene)){
    tr <- S4Vectors::split(tr, tr$entrez)
  }
  return(tr)
}

#' compute_regulatory_potential
#' 
#' @param peaks GRanges
#' @param tss GRanges or GRangesList
#' @param max_dist maximum distance between peak and tss for inclustion in regulatory 
#' potential score
#' @return vector of length tss, with regulatory potential score for each tss
#' @export
#' @author Alicia Schep
compute_regulatory_potential <- function(peaks, tss, max_dist = 100000){
  stopifnot(inherits(peaks,"GRanges"))
  score_func <- function(x){
      sum(exp(-(0.5 + 4 * x / 100000)))
  }
  out <- rep(0, length(tss))
  if (inherits(tss,"GRangesList")){
    names(tss) = 1:length(tss)
    tss_flat <- unlist(tss)
    tss_group <- as.numeric(names(tss_flat))
    score_df <- as.data.frame(GenomicRanges::findOverlaps(peaks, 
                                                   tss_flat, 
                                                   maxgap = max_dist, 
                                                   ignore.strand = TRUE)) %>%
      mutate(group = tss_group[subjectHits], 
             tmp_dists = GenomicRanges::distance(peaks[queryHits],
                                                 tss_flat[subjectHits],
                                                 ignore.strand = TRUE)) %>%
      group_by(group, queryHits) %>%
      summarise(dists = min(tmp_dists)) %>% 
      summarise(scores = score_func(dists))
    out[score_df$group] = score_df$scores
  } else if (inherits(tss, "GRanges")){
    score_df <- as.data.frame(GenomicRanges::findOverlaps(peaks, 
                                                   tss, 
                                                   maxgap = max_dist, 
                                                   ignore.strand = TRUE)) %>%
      mutate(dists = GenomicRanges::distance(peaks[queryHits],tss[subjectHits], ignore.strand = TRUE)) %>%
      group_by(subjectHits) %>% 
      summarise(scores = score_func(dists))
    out[score_df$subjectHits] = score_df$scores
  } else{
    stop("Invalid tss input, must be GRangesList or GRanges")
  }
  return(out)
}


