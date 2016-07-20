
##' plot genomic coverage along the gene
##' Adapted from gChipseq function of the same name. Credit to Jinfeng Liu
##' plot genomic coverage along the gene
##'
##' @param grl GRangesList: typically the results of importing wig/bw/... files for each sample
##' @param dat.exon data.frame: exon information
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
plot_track_view <- function(grl, dat.exon, target.range, genome, ymax, symbol,
                             bg.title='grey50', colors=NULL, sync=FALSE,
                             type = NULL, scale.group=1, hm.thresh=4,
                             scaling=NULL, hm.binsize = 1000, showSNPs=FALSE,
                             snp.gr = NULL, tss.gr = NULL,
                             dtrack.kwargs=list(), gtrack.kwargs=list(),
                             grtrack.kwargs=list(), snptrack.kwargs=list(), ...) {
  # Get the chromosome
  chr <- as.character(seqnames(target.range))
  options(ucscChromosomeNames = FALSE)
  
  if(missing(symbol)){
    main.title <- paste0("Chr", chr, ':', start(target.range), "-",
                         end(target.range))
  } else {
    main.title <- symbol
  }
  annotation.title <- "Annotation"
  
  # Decide the type of datatrack to plot if not provided
  if(is.null(type)){
    type <- ifelse(length(grl)>hm.thresh, "heatmap", "hist")
  }
  
  gtrack.defaults <- list()
  grtrack.defaults <- list(fill='lightslateblue', background.title='lightslateblue')
  snptrack.defaults <- list(fill='black', background.title='black', chromosome=chr,
                            genome=genome, cex.title=0.6)
  if(type=='hist'){
    # Make colors for the plots
    if(is.null(colors)){
      if(length(grl)>8){
        colors <- rep(RColorBrewer::brewer.pal(8, 'Dark2'), 
                      length.out=length(grl))
      } else {
        colors <- RColorBrewer::brewer.pal(length(grl), 'Dark2')
      }
    }
    colors <- rep(colors, length.out=length(grl))
    names(colors) <- names(grl)
  } else if(type=='heatmap'){
    # Shrink the genome track so the heatmap is more prevalent
    grtrack.defaults <- modifyList(grtrack.defaults, list(size=0.5))
    
  }
  
  # Override defaults with arguments
  gtrack.params <- modifyList(gtrack.defaults, gtrack.kwargs)
  grtrack.params <- modifyList(grtrack.defaults, grtrack.kwargs)
  snptrack.params <- modifyList(snptrack.defaults, snptrack.kwargs)

  # Scale data range
  if ( missing(ymax) ) {
    score.max <- sapply(grl, function(x) { max(score(x)) } )
    if ( sync ) {
      scale.group <- rep(scale.group, length.out=length(grl))
      grps <- unique(scale.group)
      for (g in grps) {
        index.g <- which(scale.group == g)
        score.max[index.g] <- max(score.max[index.g])
      }
      names(score.max) <- names(grl)
    }
  } else {
    score.max <- rep(ymax, length.out=length(grl))
    names(score.max) <- names(grl)
  }
  
    # Add title
  bg.title <- rep(bg.title, length.out=length(grl))
  names(bg.title) <- names(grl)
  
  dtrack <- make_data_tracks(grl, target.range, genome, chr = chr, bg.title, score.max = score.max, 
                           type=type, dtrack.kwargs = dtrack.kwargs, colors = colors,
                           scaling = scaling, binwidth = hm.binsize)
  
  # Add genome tracks
  gtrack <- GenomeAxisTrack()
  displayPars(gtrack) <- gtrack.params
  
  grtrack <- GeneRegionTrack(dat.exon, genome=genome, chromosome=chr, 
                             name=annotation.title, transcriptAnnotation="symbol")
  displayPars(grtrack) <- grtrack.params
  
  # Match the seqname style and remove strand info on target.range
  if(showSNPs){
    if(is.null(snp.gr)){
      stop("Please supply a SNP database object")
    }
    snptrack <- make_snp_track(snp.gr, target.range, strack.kwargs = snptrack.params)
  } else{
    snptrack=NULL
  }
  
  if(!is.null(tss.gr)){
    tssTrack <- make_tss_track(tss.gr, target.range, chr)
  } else{
    tssTrack <- NULL
  }
  
  # Add tracks
  tracklist <- list()
  tracklist <- c(gtrack, grtrack, tssTrack, dtrack, snptrack)
  
  # Confirm that all tracks have the same chromosome before plotting
  chr_list <- unlist(lapply(tracklist, function(x) 
    {if("chromosome" %in% slotNames(x)){return(x@chromosome)}}))
  if(length(unique(chr_list))!=1){
    warning("Tracks have different chromosomes. Visualization may not be complete")
  }
  
  # Plot tracks
  ptracks <- plotTracks(tracklist, main=main.title, chromosome = chr, 
                        from=start(target.range), to=end(target.range),...)
  return(invisible(ptracks))
}

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
  if (as.character(strand(gr)) == '+') {
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

bin_snp_in_range <- function(snpGR, gr, binwidth=1000){
  
  tiled.range <- make_tiled_range(gr, binwidth=binwidth)
  overlap <- findOverlaps(snpGR, tiled.range)
  bin.split <- splitAsList(snpGR,factor(subjectHits(overlap)))
  bin.score <- lapply(bin.split, function(x) length(x))
  mcols(tiled.range)[['snpCount']] <- unlist(bin.score)
  return(tiled.range)
}

make_tiled_range <- function(gr, binwidth){
  # Bin Range, currently only supports 1 range
  tiled.range <- tile(gr, width=binwidth)[[1]]
  
  # Remove strand information
  strand(tiled.range) <- "*"
  
  return(tiled.range)
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

#' make_data_tracks
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
make_data_tracks <- function(cvg.L, gr, genome, chr, bg.title, score.max, colors = NULL,
                           type=NULL, binwidth=1000, val="score", scaling=NULL,
                           dtrack.kwargs=list()){
  # Compile sample coverages as datatracks
  dtrack <- list()
  if(type == 'heatmap'){
    hm.gradient <- colorRampPalette(brewer.pal(9, "BuGn"))(100)
    binned.cvg <- bin_coverage_in_range(cvg.L, gr, binwidth, val, scaling)
    heatmap.params <- list(range=binned.cvg, type='heatmap', chromosome=chr, genome=genome,
                           background.title='gray50', showSampleNames = TRUE,
                           cex.sampleNames = 0.6, cex.axis=0.6, gradient = hm.gradient)
    heatmap.params <- modifyList(heatmap.params, dtrack.kwargs)
    dtrack[[1]] <- do.call(Gviz::DataTrack, heatmap.params)
  } else if(type == 'hist'){
    for (g in names(cvg.L)) {
      coverage.params <- list(range=cvg.L[[g]], chromosome=chr, type='hist',
                              ylim=c(0,score.max[g]),background.title=bg.title[g],
                              genome=genome, name=g, col.histogram=colors[g],
                              fill.histogram=colors[g])
      coverage.params <- modifyList(coverage.params, dtrack.kwargs)
      dtrack[[g]] <- do.call(Gviz::DataTrack, coverage.params)
    }
    names(dtrack) = names(cvg.L)
  } else {
    stop('No valid type supplied for datatrack')
  }
  return(dtrack)
}

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

make_snp_track <- function(snp.gr, target.range, strack.kwargs=list()){
  seqlevelsStyle(target.range) <- seqlevelsStyle(snp.gr)
  strand(target.range) <- "*"
  snp_hits <- findOverlaps(snp.gr, target.range, type='within')
  snp_locs <- snp.gr[queryHits(snp_hits)]
  sorted_snps <- sort_snps(snp_locs)
  grouping <- names(mcols(sorted_snps))
  strack.params <- list(range = sorted_snps, name='snps', showAxis=FALSE,
                        groups = grouping, type='p', legend=FALSE, rot.title=0,
                        size = 2)
  strack.params <- modifyList(strack.params, strack.kwargs)
  strack <- do.call(Gviz::DataTrack, strack.params)
  return(strack)
}

make_tss_track <- function(tss, target.range, chr){
  tss_overlap <- findOverlaps(tss, target.range)
  if(length(tss_overlap)==0){
    return(NULL)
  }
  tss_in_range <- tss[queryHits(tss_overlap)]
  tssTrack <- AnnotationTrack(tss_in_range, shape='arrow', name = 'TSS',
                              id = mcols(tss_in_range)['names'], chromosome=chr,
                              showFeatureId=TRUE) 
  return(tssTrack)
  }