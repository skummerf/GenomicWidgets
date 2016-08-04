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
  tiled_range <- tile(gr, width=binwidth)[[1]]
  
  # Remove strand information
  strand(tiled_range) <- "*"
  
  return(tiled_range)
}

# =============================================================================
# TSS track and helper functions

#' Make TSS Track
#'
#' @param tss 
#' @param target_range 
#' @param chr 
#' @author Justin Finkle
#' @return
#' @export
#'
#' @examples
make_tss_track <- function(tss, target_range, chr, kwargs=list()){
  # Identify TSS in the target region
  tss_overlap <- findOverlaps(tss, target_range)
  if(length(tss_overlap)==0){
    return(NULL)
  }
  tss_in_range <- tss[queryHits(tss_overlap)]
  
  # TSS are all ranges of width 1, so to make the arrows have a size the ranges
  # are adjusted
  tss_df <- as.data.frame(tss_in_range)
  tss_df[['start']] <- tss_df[['start']] - round(width(target_range)*0.05)*(tss_df[['strand']]=='-')
  tss_df[['end']] <- tss_df[['end']] + round(width(target_range)*0.05)*(tss_df[['strand']]=='+')
  tss_in_range <- GRanges(tss_df)
  
  tss_defaults <- list(range=tss_in_range, 
                       shape='arrow', 
                       name = 'TSS',
                       id = mcols(tss_in_range)[['names']], 
                       chromosome=chr,
                       showFeatureId=TRUE, 
                       cex=0.75, 
                       fontcolor.feature='black',
                       rot.title = 0, 
                       background.title='black') 
  tss_params <- modifyList(tss_defaults, kwargs)
  tssTrack <- do.call(Gviz::AnnotationTrack, tss_params)
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
  grtrack_defaults <- list(range = exon_data, 
                           fill = 'lightslateblue', 
                           rot.title = 0,
                           cex.title = 0.75,
                           background.title = 'lightslateblue', 
                           name = 'Transcripts', 
                           transcriptAnnotation = "transcript",
                           chromosome=chr, 
                           genome=genome)
  grtrack_params <- modifyList(grtrack_defaults, kwargs)
  grtrack <- do.call(Gviz::GeneRegionTrack, grtrack_params)
  # Preemptively get stacks to see if the annotations should be collapsed
  grtrack <- Gviz:::setStacks(grtrack)
  if(max(grtrack@stacks) > 3 & is.null(kwargs$stacking)){
    grtrack_params <- modifyList(grtrack_params, list(stacking='dense'))
    grtrack <- do.call(Gviz::GeneRegionTrack, grtrack_params)
  }
  return(grtrack)
}
# =============================================================================
# GenomeAxis Track and helper functions

#' Make GenomeAxisTrack
#'
#' @param target_range 
#' @param kwargs 
#'
#' @return
#' @export
#'
#' @examples
make_gatrack <- function(kwargs = list()){
  gatrack_defaults <- list()
  gatrack_params <- modifyList(gatrack_defaults, kwargs)
  gatrack <- do.call(Gviz::GenomeAxisTrack, gatrack_params)
  return(gatrack)
}


# =============================================================================
# SNP track and helper functions

#' Make SNP Track
#'
#' @param snp_gr 
#' @param target_range 
#' @param strack.kwargs 
#' @author Justin Finkle
#' @return
#' @export
#'
#' @examples
make_snp_track <- function(snp_gr, target_range, kwargs=list()){
  seqlevelsStyle(target_range) <- seqlevelsStyle(snp_gr)
  strand(target_range) <- "*"
  snp_hits <- findOverlaps(snp_gr, target_range, type='within')
  snp_locs <- snp_gr[queryHits(snp_hits)]
  if(length(snp_locs)==0){
    return(NULL)
  }
  sorted_snps <- sort_snps(snp_locs)
  grouping <- names(mcols(sorted_snps))
  strack_defaults <- list(range = sorted_snps, 
                          name='SNPs', 
                          showAxis=FALSE,
                          groups = grouping, 
                          type='p', 
                          legend=FALSE,
                          background.title='black',
                          rot.title = 0,
                          showSampleNames=FALSE,
                          cex.sampleNames = 0.6)
  strack_params <- modifyList(strack_defaults, kwargs)
  strack <- do.call(Gviz::DataTrack, strack_params)
  return(strack)
}

#' Title
#'
#' @param snp_gr 
#' @param type_col 
#' @author Justin Finkle
#' @return
#' @export
#'
#' @examples
sort_snps <- function(snp_gr, type_col = 'CONTEXT'){
  snp_types <- mcols(snp_gr)[[type_col]]
  for(type in unique(snp_types)){
    mcols(snp_gr)[type] <- (mcols(snp_gr)[[type_col]]==type)*which(type==unique(snp_types))
  }
  new_df <- as.data.frame(mcols(snp_gr)[, names(mcols(snp_gr)) 
                                        %in% unique(snp_types)])
  names(new_df) <- gsub("_variant", "", names(new_df))
  new_df[new_df == 0 ] <- NA
  mcols(snp_gr) <- new_df
  return(snp_gr)
}

# =============================================================================
# Data Tracks and helper functions

#' Make Data Tracks
#' @description Builds Gviz DataTrack objects for plotting based on coverage
#'
#' @param gr GRange to plot
#' @param genome str
#' @param chr str
#' @param bg_title str: color for track panel backgrounds 
#' @param colors list-optional: colors for coverage tracks
#' @param type str: 'hist' for multiple coverage tracks, 'heatmap' for condensed view
#' @param binwidth scalar
#' @param val str: value in cvg.L elements to use for binning coverage. Default is "score"
#' @param cvg_list 
#' @param kwargs 
#' @param scaling 
#'
#' @return dtrack list: DataTrack objects
#' @export
#' @author Justin Finkle
#' @import Gviz
#'
#' @examples
make_data_tracks <- function(cvg_list, gr, genome, chr,
                             bg_title = 'black', 
                             colors = NULL, 
                             type=NULL, 
                             binwidth=1000, 
                             val="score", 
                             hm_scaling = log,
                             kwargs=list()){
  # Compile sample coverages as datatracks
  dtrack <- list()
  if(type == 'heatmap'){
    hm_gradient <- colorRampPalette(brewer.pal(9, "BuGn"))(100)
    binned_cvg <- bin_coverage_in_range(cvg_list = cvg_list, gr = gr,
                                        binwidth = binwidth,
                                        val = val, scaling = hm_scaling)
    hm_defaults <- list(range = binned_cvg, 
                        type = 'heatmap', 
                        chromosome = chr,
                        genome = genome, 
                        background.title = bg_title, 
                        cex.axis = 0.6,
                        showSampleNames = TRUE, 
                        cex.sampleNames = 0.6,
                        gradient = hm_gradient,
                        showTitle = FALSE)
    hm_params <- modifyList(hm_defaults, kwargs)
    dtrack[[1]] <- do.call(Gviz::DataTrack, hm_params)
  } else if(type == 'hist'){
    
    # Make colors for the plots
    if(is.null(colors)){
      if(length(cvg_list)>8){
        colors <- rep(RColorBrewer::brewer.pal(8, 'Dark2'), 
                      length.out=length(cvg_list))
        } else {
          colors <- RColorBrewer::brewer.pal(length(cvg_list), 'Dark2')
        }
      }
    colors <- rep(colors, length.out=length(cvg_list))
    names(colors) <- names(cvg_list)
    
    # Add title
    bg_title <- rep(bg_title, length.out=length(cvg_list))
    names(bg_title) <- names(cvg_list)
    
    # Scale the data for the histograms
    score_max <- max(sapply(cvg_list, function(x) { max(score(x)) } ))
    
    for (g in names(cvg_list)) {
      coverage_defaults <- list(range=cvg_list[[g]], 
                                chromosome=chr, 
                                type='hist',
                                ylim=c(0, score_max),
                                background.title=bg_title[g],
                                genome=genome, name=g, 
                                col.histogram=colors[g],
                                fill.histogram=colors[g])
      coverage_params <- modifyList(coverage_defaults, kwargs)
      dtrack[[g]] <- do.call(Gviz::DataTrack, coverage_params)
    }
    names(dtrack) = names(cvg_list)
  } else {
    stop('No valid type supplied for datatrack')
  }
  return(dtrack)
}

bin_coverage_in_range <- function(cvg_list, gr, 
                                  binwidth=1000, 
                                  val="score", 
                                  weighted = TRUE,
                                  scaling = NULL){
  
  if(is.null(scaling)){
    scaling <- function(x){x}
  }
  
  # Bin Range, currently only supports 1 range
  tiled_range <- tile(gr, width=binwidth)[[1]]
  
  # Remove strand information
  strand(tiled_range) <- "*"
  
  for(g in names(cvg_list)){
    # Find overlap between the coverage and the tiled range
    overlap <- findOverlaps(cvg_list[[g]], tiled_range)
    if(weighted){
      wscores <- width(cvg_list[[g]])*mcols(cvg_list[[g]])[[val]]
    } else{
      wscores <- mcols(cvg_list[[g]])[[val]]
    }
    # Split scores from coverage range into their appropriate bin and sum
    bin_split <- splitAsList(wscores, factor(subjectHits(overlap)))
    
    # Need a better way to scale the data
    bin_score <- lapply(bin_split, function(x) sum(x))
    mcols(tiled_range)[[g]] <- unlist(bin_score)
  }
  return(tiled_range)
}

# =============================================================================
# Deprecated Functions

# This is based on code from plotGeneCoverage in gChipseq. It seems unnecessary
# to provide this much flexibility in scaling
scale_track_data <- function(cvg_list, ymax, sync=FALSE, scale_group=1){
  # Scale data range
  if (missing(ymax)) {
    score_max <- sapply(cvg_list, function(x) { max(score(x)) } )
    if (sync) {
      scale_group <- rep(scale_group, length.out=length(cvg_list))
      grps <- unique(scale_group)
      for (g in grps) {
        index.g <- which(scale_group == g)
        score_max[index.g] <- max(score_max[index.g])
      }
      names(score_max) <- names(cvg_list)
    }
  } else {
    score_max <- rep(ymax, length.out=length(cvg_list))
    names(score_max) <- names(cvg_list)
  }
  return(score_max)
}

# Only used if want to display snps as a heatmap
bin_snp_in_range <- function(snp_gr, gr, binwidth=1000){
  tiled_range <- make_tiled_range(gr, binwidth=binwidth)
  overlap <- findOverlaps(snp_gr, tiled_range)
  bin_split <- splitAsList(snp_gr,factor(subjectHits(overlap)))
  bin_score <- lapply(bin_split, function(x) length(x))
  mcols(tiled_range)[['snpCount']] <- unlist(bin_score)
  return(tiled_range)
}

