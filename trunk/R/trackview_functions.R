
##' plot genomic coverage along the gene
##' Adapted from gChipseq function of the same name. Credit to Jinfeng Liu
##' plot genomic coverage along the gene
##'
##' @param cvg.L GRangesList: typically the results of importing wig/bw/... files for each sample
##' @param transcripts.gr data.frame: exon information
##' @param target.range GRange: range to plot 
##' @param genome str: genome build, e.g. hg19, GRCh38, etc.
##' @param ymax vector-optional: ymax for the data tracks
##' @param symbol str-optional: symbol to plot
##' @param bg.title: str-optional: background color for the title panel
##' @param colors str-optional: colors for coverage data tracks
##' @param sync logical-optional: whether to sync the ymax for all data tracks
##' @param scale.group scalar-optional: groups to use the same y scale
##' @param hm.thresh integer-optional: threshold above which a heatmap is plotted. Can change threshold or override with type
##' @param type str-optional: type of datatrack to plot. Currently only 'hist' and 'heatmap' are supported
##' @param dtrack.kwargs list-optional: additional arguments for DataTrack. These should only be used for global arguments that will apply to all datatracks
##' @param gtrack.kwargs list-optional: additional arguments for GenomeAxisTrack
##' @param grtrack.kwargs list-optional: additional arguments for GeneRegionTrack
##' @param snptrack.kwargs list-optional: additional arguments for the SNP DataTrack
##'
##' @return ptracks list: track objects plotted by Gviz
##' @import Gviz gChipseq RColorBrewer
##' @export
##' @author Justin Finkle
plot_track_view <- function(cvg.L, target.range, transcripts.gr, genome, ymax, symbol,
                            bg.title = 'black', colors = NULL, sync = FALSE,
                            type = NULL, scale.group = 1, hm.thresh = 4,
                            scaling = NULL, hm.binsize = 1000, snp.gr = NULL, 
                            tss.gr = NULL, ucsc = FALSE, gatrack.kwargs=list(),
                            grtrack.kwargs=list(), tsstrack.kwargs=list(), 
                            snptrack.kwargs=list(), dtrack.kwargs=list(), 
                            plot.kwargs = list()) {
  # Initialize
  chr <- as.character(seqnames(target.range))
  options(ucscChromosomeNames = ucsc)
  main.title <- make_title(target.range, chr)
  snptrack <- NULL
  tssTrack <- NULL
  
  # Decide the type of datatrack to plot if not provided
  if(is.null(type)){
    type <- ifelse(length(cvg.L)>hm.thresh, "heatmap", "hist")
  }
  
  # Make tracks
  dtrack <- make_data_tracks(cvg.L, target.range, genome, chr = chr, bg.title, 
                             ymax = ymax, scale.group = scale.group, sync = sync,
                             type=type, kwargs = dtrack.kwargs, colors = colors,
                             scaling = scaling, binwidth = hm.binsize)
  
  # Add genome tracks
  gtrack <- make_gatrack(kwargs = gatrack.kwargs)
  grtrack <- make_grtrack(transcripts.gr, genome, chr, kwargs = grtrack.kwargs)
  
  if(!is.null(snp.gr)){
    snptrack <- make_snp_track(snp.gr, target.range, kwargs = snptrack.kwargs)
  }
  
  if(!is.null(tss.gr)){
    tssTrack <- make_tss_track(tss.gr, target.range, chr, kwargs = tsstrack.kwargs)
  }
  
  # Add tracks
  tracklist <- list()
  tracklist <- c(gtrack, grtrack, tssTrack, snptrack, dtrack)
  
  # Confirm that all tracks have the same chromosome before plotting
  chr_list <- unlist(lapply(tracklist, chr_check))
  if(length(unique(chr_list))!=1){
    warning("Tracks have different chromosomes. Visualization may not be complete")
  }
  
  # Plot tracks
  plot.defaults <- list(trackList = tracklist, main = main.title, chromosome = chr,
                        from = start(target.range), to = end(target.range),
                        cex.main = 1)
  plot.params <- modifyList(plot.defaults, plot.kwargs)
  ptracks <- do.call(Gviz::plotTracks, plot.params)
  return(invisible(ptracks))
}

make_title <- function(target.range, chr){
  title <- paste0("Chr", chr, ':', start(target.range), "-", end(target.range))
  return(title)
}

chr_check <- function(track){
  if("chromosome" %in% slotNames(track)){
    return(track@chromosome)
  }
}