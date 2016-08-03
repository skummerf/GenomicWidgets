#' Title
#'
#' @param target_range 
#' @param cvg_files 
#' @param genome 
#' @param db_object 
#' @param tx_data 
#' @param bg_title 
#' @param colors 
#' @param type 
#' @param hm_thresh 
#' @param cvg_scaling 
#' @param hm_binsize 
#' @param snp_gr 
#' @param tss_gr 
#' @param ucsc 
#' @param gatrack_kwargs 
#' @param grtrack_kwargs 
#' @param tsstrack_kwargs 
#' @param snptrack_kwargs 
#' @param dtrack_kwargs 
#' @param plot_kwargs 
#'
#' @return
#' @export
#'
#' @examples
make_track_function <- function(target_range, cvg_files, genome, db_object, tx_data,
                                cvg_scaling = NULL,
                                bg_title = 'black', 
                                colors = NULL, 
                                type = NULL, 
                                hm_thresh = 4,
                                hm_binsize = 1000,
                                hm_scaling = NULL,
                                snp_gr = NULL, 
                                tss_gr = NULL, 
                                ucsc = FALSE, 
                                gatrack_kwargs=list(),
                                grtrack_kwargs=list(), 
                                tsstrack_kwargs=list(), 
                                snptrack_kwargs=list(), 
                                dtrack_kwargs=list(), 
                                plot_kwargs = list()){
  plot_tv <- function(range){
    exon_data <- get_tx_annotation(db_object = txdb, 
                                   range = view_range, 
                                   tx_data = tx_data,
                                   no_introns=TRUE)
    cvg_list <- get_coverage_in_range(bwList = cvg_files,
                                      target_range = target_range, 
                                      names = names(cvg_files), 
                                      cvg_scaling = cvg_scaling)
    tl <- plot_track_view(cvg_list, target_range, exon_data, genome,
                          bg_title = bg_title,
                          colors = colors, 
                          type = type, 
                          hm_thresh = hm_thresh,
                          hm_binsize = hm_binsize,
                          hm_scaling = hm_scaling,
                          snp_gr = snp_gr, 
                          tss_gr = tss_gr, 
                          ucsc = ucsc, 
                          gatrack_kwargs = gatrack_kwargs,
                          grtrack_kwargs = grtrack_kwargs,
                          tsstrack_kwargs = tsstrack_kwargs,
                          snptrack_kwargs = snptrack_kwargs, 
                          dtrack_kwargs = dtrack_kwargs, 
                          plot_kwargs = plot_kwargs)
  }
  return(plot_tv)
}


##' plot genomic coverage along the gene
##' Adapted from gChipseq function of the same name. Credit to Jinfeng Liu
##' plot genomic coverage along the gene
##'
##' @param cvg_list GRangesList: typically the results of importing wig/bw/... files for each sample
##' @param target_range GRange: range to plot 
##' @param exon_data data.frame: exon information
##' @param genome str: genome build, e.g. hg19, GRCh38, etc.
##' @param bg_title 
##' @param colors str-optional: colors for coverage data tracks
##' @param type str-optional: type of datatrack to plot. Currently only 'hist' and 'heatmap' are supported
##' @param hm_thresh integer-optional: threshold above which a heatmap is plotted. Can change threshold or override with type
##' @param scaling 
##' @param hm_binsize 
##' @param snp_gr 
##' @param tss_gr 
##' @param ucsc 
##' @param gatrack_kwargs 
##' @param grtrack_kwargs list-optional: additional arguments for GeneRegionTrack
##' @param tsstrack_kwargs 
##' @param snptrack_kwargs list-optional: additional arguments for the SNP DataTrack
##' @param dtrack_kwargs list-optional: additional arguments for DataTrack. These should only be used for global arguments that will apply to all datatracks
##' @param plot_kwargs 
##'
##' @return ptracks list: track objects plotted by Gviz
##' @import Gviz gChipseq RColorBrewer
##' @export
##' @author Justin Finkle
plot_track_view <- function(cvg_list, target_range, exon_data, genome,
                            bg_title = 'black', 
                            colors = NULL, 
                            type = NULL, 
                            hm_thresh = 4,
                            hm_binsize = 1000, 
                            hm_scaling = NULL,
                            snp_gr = NULL, 
                            tss_gr = NULL, 
                            ucsc = FALSE, 
                            gatrack_kwargs=list(),
                            grtrack_kwargs=list(), 
                            tsstrack_kwargs=list(), 
                            snptrack_kwargs=list(), 
                            dtrack_kwargs=list(), 
                            plot_kwargs = list()) {
  # Initialize
  chr <- as.character(seqnames(target_range))
  options(ucscChromosomeNames = ucsc)
  main_title <- make_title(target_range = target_range, 
                           chr = chr)
  snptrack <- NULL
  tssTrack <- NULL
  
  # Decide the type of datatrack to plot if not provided
  if(is.null(type)){
    type <- ifelse(length(cvg_list) > hm_thresh, "heatmap", "hist")
  }
  
  # Make tracks
  dtrack <- make_data_tracks(cvg_list = cvg_list, 
                             gr = target_range, 
                             genome = genome, 
                             chr = chr, 
                             bg_title = bg_title, 
                             type=type, 
                             colors = colors,
                             binwidth = hm_binsize,
                             hm_scaling = hm_scaling,
                             kwargs = dtrack_kwargs)
  
  # Add genome tracks
  gatrack <- make_gatrack(kwargs = gatrack_kwargs)
  grtrack <- make_grtrack(exon_data = exon_data, 
                          genome = genome, 
                          chr = chr, 
                          kwargs = grtrack_kwargs)
  
  if(!is.null(snp_gr)){
    snptrack <- make_snp_track(snp_gr = snp_gr,
                               target_range = target_range, 
                               kwargs = snptrack_kwargs)
  }
  
  if(!is.null(tss_gr)){
    tssTrack <- make_tss_track(tss = tss_gr,
                               target_range = target_range,
                               chr = chr, 
                               kwargs = tsstrack_kwargs)
  }
  
  # Add tracks
  tracklist <- list()
  tracklist <- c(GeneAxis=gatrack, GeneRegion=grtrack, 
                 TSS=tssTrack, SNP=snptrack, dtrack)
  
  # Confirm that all tracks have the same chromosome before plotting
  chr_list <- unlist(lapply(tracklist, chr_check))
  if(length(unique(chr_list))!=1){
    warning("Tracks have different chromosomes. Visualization may not be complete")
  }
  
  # Plot tracks
  plot_defaults <- list(trackList = tracklist, 
                        main = main_title, 
                        chromosome = chr,
                        from = start(target_range), 
                        to = end(target_range),
                        cex.main = 1)
  plot.params <- modifyList(plot_defaults, plot_kwargs)
  ptracks <- do.call(Gviz::plotTracks, plot.params)
  return(invisible(ptracks))
}

# ========================= Minor Helper Functions =============================
make_title <- function(target_range, chr){
  title <- paste0("Chr", chr, ':', start(target_range), "-", end(target_range))
  return(title)
}

chr_check <- function(track){
  if("chromosome" %in% slotNames(track)){
    return(track@chromosome)
  }
}