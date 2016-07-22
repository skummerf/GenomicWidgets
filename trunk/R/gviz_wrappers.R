# =============================================================================
# =============================================================================
# These functions are intended to provide customizable wrappers for the GViz
# tracks. Logical defaults for use in chipVis are provided, but each function
# can also take a list of additional arguments that are passed to the core GViz
# function
# =============================================================================
# =============================================================================


# =============================================================================
# General Functions
make_tiled_range <- function(gr, binwidth){
  # Bin Range, currently only supports 1 range
  tiled.range <- tile(gr, width=binwidth)[[1]]
  
  # Remove strand information
  strand(tiled.range) <- "*"
  
  return(tiled.range)
}

# =============================================================================
# TSS track and helper functions

#' Make TSS Track
#'
#' @param tss 
#' @param target.range 
#' @param chr 
#' @author Justin Finkle
#' @return
#' @export
#'
#' @examples
make_tss_track <- function(tss, target.range, chr, kwargs=list()){
  # Identify TSS in the target region
  tss_overlap <- findOverlaps(tss, target.range)
  if(length(tss_overlap)==0){
    return(NULL)
  }
  tss_in_range <- tss[queryHits(tss_overlap)]
  
  # TSS are all ranges of width 1, so to make the arrows have a size the ranges
  # are adjusted
  tss_df <- as.data.frame(tss_in_range)
  tss_df[['start']] <- tss_df[['start']] - round(width(target.range)*0.05)*(tss_df[['strand']]=='-')
  tss_df[['end']] <- tss_df[['end']] + round(width(target.range)*0.05)*(tss_df[['strand']]=='+')
  tss_in_range <- GRanges(tss_df)
  
  tss.defaults <- list(range=tss_in_range, shape='arrow', name = 'TSS',
                       id = mcols(tss_in_range)[['names']], chromosome=chr,
                       showFeatureId=TRUE, cex=0.75, fontcolor.feature='black',
                       rot.title = 0, background.title='black') 
  tss.params <- modifyList(tss.defaults, kwargs)
  tssTrack <- do.call(Gviz::AnnotationTrack, tss.params)
  return(tssTrack)
}
# =============================================================================
# GeneRegionTrack and helper functions

#' Make GeneRegionTrack
#'
#' @param exon_data 
#' @param genome 
#' @param chr 
#' @param kwargs 
#'
#' @return grtrack
#' @export
#'
#' @examples
make_grtrack <- function(exon_data, genome, chr, kwargs = list()){
  grtrack.defaults <- list(range = exon_data, fill = 'lightslateblue', 
                           rot.title = 0, background.title = 'lightslateblue', 
                           name = 'Transcripts', transcriptAnnotation = "symbol",
                           chromosome=chr, genome=genome)
  grtrack.params <- modifyList(grtrack.defaults, kwargs)
  grtrack <- do.call(Gviz::GeneRegionTrack, grtrack.params)
  # Preemptively get stacks to see if the annotations should be collapsed
  grtrack <- Gviz:::setStacks(grtrack)
  if(max(grtrack@stacks) > 3 & is.null(kwargs$stacking)){
    grtrack.params <- modifyList(grtrack.params, list(stacking='dense'))
    grtrack <- do.call(Gviz::GeneRegionTrack, grtrack.params)
  }
  return(grtrack)
}
# =============================================================================
# GenomeAxis Track and helper functions

#' Make GenomeAxisTrack
#'
#' @param target.range 
#' @param kwargs 
#'
#' @return
#' @export
#'
#' @examples
make_gatrack <- function(kwargs = list()){
  gatrack.defaults <- list()
  gatrack.params <- modifyList(gatrack.defaults, kwargs)
  gatrack <- do.call(Gviz::GenomeAxisTrack, gatrack.params)
  return(gatrack)
}


# =============================================================================
# SNP track and helper functions

#' Make SNP Track
#'
#' @param snp.gr 
#' @param target.range 
#' @param strack.kwargs 
#' @author Justin Finkle
#' @return
#' @export
#'
#' @examples
make_snp_track <- function(snp.gr, target.range, kwargs=list()){
  seqlevelsStyle(target.range) <- seqlevelsStyle(snp.gr)
  strand(target.range) <- "*"
  snp_hits <- findOverlaps(snp.gr, target.range, type='within')
  snp_locs <- snp.gr[queryHits(snp_hits)]
  if(length(snp_locs)==0){
    return(NULL)
  }
  sorted_snps <- sort_snps(snp_locs)
  grouping <- names(mcols(sorted_snps))
  strack.params <- list(range = sorted_snps, name='SNPs', showAxis=FALSE,
                        groups = grouping, type='p', legend=FALSE,
                        background.title='black', showSampleNames=TRUE,
                        cex.sampleNames = 0.6)
  strack.params <- modifyList(strack.params, kwargs)
  strack <- do.call(Gviz::DataTrack, strack.params)
  return(strack)
}

#' Title
#'
#' @param snp.gr 
#' @param type.col 
#' @author Justin Finkle
#' @return
#' @export
#'
#' @examples
sort_snps <- function(snp.gr, type.col = 'CONTEXT'){
  snp_types <- mcols(snp.gr)[[type.col]]
  for(type in unique(snp_types)){
    mcols(snp.gr)[type] <- (mcols(snp.gr)[[type.col]]==type)*which(type==unique(snp_types))
  }
  new_df <- as.data.frame(mcols(snp.gr)[, names(mcols(snp.gr)) 
                                        %in% unique(snp_types)])
  names(new_df) <- gsub("_variant", "", names(new_df))
  new_df[new_df == 0 ] <- NA
  mcols(snp.gr) <- new_df
  return(snp.gr)
}

bin_snp_in_range <- function(snpGR, gr, binwidth=1000){
  
  tiled.range <- make_tiled_range(gr, binwidth=binwidth)
  overlap <- findOverlaps(snpGR, tiled.range)
  bin.split <- splitAsList(snpGR,factor(subjectHits(overlap)))
  bin.score <- lapply(bin.split, function(x) length(x))
  mcols(tiled.range)[['snpCount']] <- unlist(bin.score)
  return(tiled.range)
}
# =============================================================================
# Data Tracks and helper functions

#' Make Data Tracks
#' @description Builds Gviz DataTrack objects for plotting based on coverage
#'
#' @param cvg.L GRangesList containing coverage score
#' @param gr GRange to plot
#' @param genome str
#' @param chr str
#' @param bg.title str: color for track panel backgrounds 
#' @param score.max scalar: scaling
#' @param colors list-optional: colors for coverage tracks
#' @param type str: 'hist' for multiple coverage tracks, 'heatmap' for condensed view
#' @param binwidth scalar
#' @param val str: value in cvg.L elements to use for binning coverage. Default is "score"
#' @param scaling 
#' @param dtrack.kwargs list-optional: additional arguments for DataTrack. These should only be used for global arguments that will apply to all datatracks
#'
#' @return dtrack list: DataTrack objects
#' @export
#' @author Justin Finkle
#' @import Gviz
#'
#' @examples
make_data_tracks <- function(cvg.L, gr, genome, chr, ymax, sync=FALSE, 
                             bg.title = 'black', scale.group=1, colors = NULL, 
                             type=NULL, binwidth=1000, val="score", scaling=NULL,
                             kwargs=list()){
  # Compile sample coverages as datatracks
  dtrack <- list()
  if(type == 'heatmap'){
    hm.gradient <- colorRampPalette(brewer.pal(9, "BuGn"))(100)
    binned.cvg <- bin_coverage_in_range(cvg.L, gr, binwidth, val, scaling)
    hm.defaults <- list(range = binned.cvg, type = 'heatmap', chromosome = chr,
                        genome = genome, background.title = bg.title, cex.axis=0.6,
                        showSampleNames = TRUE, cex.sampleNames = 0.6,
                        gradient = hm.gradient)
    hm.params <- modifyList(hm.defaults, kwargs)
    dtrack[[1]] <- do.call(Gviz::DataTrack, hm.params)
  } else if(type == 'hist'){
    
    # Make colors for the plots
    if(is.null(colors)){
      if(length(cvg.L)>8){
        colors <- rep(RColorBrewer::brewer.pal(8, 'Dark2'), 
                      length.out=length(cvg.L))
        } else {
          colors <- RColorBrewer::brewer.pal(length(cvg.L), 'Dark2')
        }
      }
    colors <- rep(colors, length.out=length(cvg.L))
    names(colors) <- names(cvg.L)
    
    # Add title
    bg.title <- rep(bg.title, length.out=length(cvg.L))
    names(bg.title) <- names(cvg.L)
    
    # Get a scaling factor for histograms
    score.max <- scale_track_data(cvg.L, ymax, sync, scale.group)
    
    for (g in names(cvg.L)) {
      coverage.defaults <- list(range=cvg.L[[g]], chromosome=chr, type='hist',
                              ylim=c(0,score.max[g]),background.title=bg.title[g],
                              genome=genome, name=g, col.histogram=colors[g],
                              fill.histogram=colors[g])
      coverage.params <- modifyList(coverage.defaults, kwargs)
      dtrack[[g]] <- do.call(Gviz::DataTrack, coverage.params)
    }
    names(dtrack) = names(cvg.L)
  } else {
    stop('No valid type supplied for datatrack')
  }
  return(dtrack)
}

bin_coverage_in_range <- function(cvg.L, gr, binwidth=1000, val="score", 
                                  scaling=NULL, weighted=TRUE){
  if(is.null(scaling)){
    scaling <- function(x){x}
  }
  # Bin Range, currently only supports 1 range
  tiled.range <- tile(gr, width=binwidth)[[1]]
  
  # Remove strand information
  strand(tiled.range) <- "*"
  
  for(g in names(cvg.L)){
    # Find overlap between the coverage and the tiled range
    overlap <- findOverlaps(cvg.L[[g]], tiled.range)
    if(weighted){
      wscores <- width(cvg.L[[g]])*mcols(cvg.L[[g]])[[val]]
    } else{
      wscores <- mcols(cvg.L[[g]])[[val]]
    }
    # Split scores from coverage range into their appropriate bin and sum
    bin.split <- splitAsList(wscores, factor(subjectHits(overlap)))
    
    # Need a better way to scale the data
    bin.score <- lapply(bin.split, function(x) scaling(sum(x)))
    elementMetadata(tiled.range)[[g]] <- unlist(bin.score)
  }
  return(tiled.range)
}

scale_track_data <- function(cvg.L, ymax, sync=FALSE, scale.group=1){
  # Scale data range
  if (missing(ymax)) {
    score.max <- sapply(cvg.L, function(x) { max(score(x)) } )
    if (sync) {
      scale.group <- rep(scale.group, length.out=length(cvg.L))
      grps <- unique(scale.group)
      for (g in grps) {
        index.g <- which(scale.group == g)
        score.max[index.g] <- max(score.max[index.g])
      }
      names(score.max) <- names(cvg.L)
    }
  } else {
    score.max <- rep(ymax, length.out=length(cvg.L))
    names(score.max) <- names(cvg.L)
  }
  return(score.max)
}
# =============================================================================