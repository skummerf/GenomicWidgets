
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
plotGeneCoverage <- function(grl, dat.exon, target.range, genome, ymax, symbol,
                             bg.title='grey50', colors=NULL, sync=FALSE, 
                             type = NULL, scale.group=1, hm.thresh=4, 
                             showSNPs=FALSE, snpDB = SNPlocs.Hsapiens.dbSNP141.GRCh38,
                             dtrack.kwargs=list(), gtrack.kwargs=list(), 
                             grtrack.kwargs=list(), snptrack.kwargs=list()) {
  # Get the chromosome
  chr <- as.character(seqnames(target.range))
  if(missing(symbol)){
    main.title <- paste0("Chr ", chr, '>', start(target.range), ":",
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
  grtrack.defaults <- list(fill='lightslateblue', cex.title=0.6,
                           background.title='lightslateblue')
  snptrack.defaults <- list(fill='black', background.title='black', type='hist',
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
    grtrack.defaults <- modifyList(grtrack.defaults, list(size=0.25))
    
  }
  
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
  
  dtrack <- makeDataTracks(grl, target.range, genome, chr, bg.title, score.max = score.max, 
                           type=type, dtrack.kwargs = dtrack.kwargs, colors = colors)
  
  # Add genome tracks
  options(ucscChromosomeNames = FALSE)
  gtrack <- GenomeAxisTrack()
  displayPars(gtrack) <- gtrack.params
  
  grtrack <- GeneRegionTrack(dat.exon, genome=genome, chromosome=chr, 
                             name=annotation.title, transcriptAnnotation="symbol")
  displayPars(grtrack) <- grtrack.params
  
  # Match the seqname style and remove strand info on target.range
  if(showSNPs){
    seqlevelsStyle(target.range) <- seqlevelsStyle(snpDB)
    strand(target.range) <- "*"
    snp_locs <- snpsByOverlaps(snpDB, target.range, type='within')
    snp_counts <- binSNPInRange(snp_locs, target.range)
    seqlevelsStyle(snp_counts) <- seqlevelsStyle(grl)
    snptrack <- DataTrack(snp_counts, name ='# SNPs')
    displayPars(snptrack) <- snptrack.params
  } else{
    snptrack=NULL
  }
    # Add tracks
  tracklist <- list()
  tracklist <- c(gtrack,grtrack, dtrack, snptrack)
  
  # Plot tracks
  ptracks <- plotTracks(tracklist, main=main.title, chromosome = chr)
  return(ptracks)
}

#' Title
#'
#' @param dat.exon 
#'
#' @author Justin Finkle
#' @return
#' @export
#'
#' @examples
getGeneRange <- function(dat.exon){
  # Make the range
  range.gene <- GRanges(seqnames=dat.exon[1,"chr"],
                        ranges=IRanges(start=min(dat.exon$start),
                                       end=max(dat.exon$end)),
                        strand=dat.exon[1,"strand"])
  return(range.gene)
}

#' Title
#'
#' @param gr
#' @param extend 
#'
#' @return
#' @export
#' @author Justin Finkle
#'
#' @examples
extendGRange <- function(gr, extend=0){
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


#' Title
#'
#' @param bwList 
#' @param target.range 
#' @param names 
#'
#' @return
#' @export
#'
#' @examples
getCoverageInRange <- function(bwList, target.range, names, scaling.factor=NULL){
  # Import only the range that matches target.range
  cvg <- lapply(bwList, function(x) import.bw(x, which=target.range))
  cvg <- GRangesList(cvg)
  if(!is.null(scaling.factor)){
    for(g in 1:length(cvg)){
      mcols(cvg[[g]])[['score']] <- mcols(cvg[[g]])[['score']]/scaling.factor[[g]]
    }
  }
  
  
  # Name the list
  if (!missing(names)){
    names(cvg) <- names
  }
  return(cvg)
}

##' get exon information from IGIS
##'
##' get exon information from IGIS
##' @param ids vector of ids to retrive the exon information
##' @param src source of id, gene|symbol
##' @param org organism, human or mouse
##' @param ... additional argument passed to igisGetmart()
##' @return a data frame of exon information
##' @import biomaRt gChipseq
##' @export
##' @author Jinfeng Liu
igisExonlist <- function(target, src=c('symbol','gene'), genome,
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
getViewRange <- function(symbol, org, genome, chr=NULL, start=NULL, end=NULL, 
                         strand=NULL){
  if(!missing(symbol)){
    exon.df <- chipVis::igisExonlist(target=symbol, org=org, genome=genome)
    view.range <- getGeneRange(exon.df)
  } else {
    view.range <- GRanges(Rle(c(chr)), IRanges(start=start, end=end), 
                          strand = strand)
  }
  return(view.range)
}

binSNPInRange <- function(snpGR, gr, binwidth=1000){
  
  tiled.range <- makeTiledRange(gr, binwidth=binwidth)
  overlap <- findOverlaps(snpGR, tiled.range)
  bin.split <- splitAsList(snpGR,factor(subjectHits(overlap)))
  bin.score <- lapply(bin.split, function(x) length(x))
  mcols(tiled.range)[['snpCount']] <- unlist(bin.score)
  return(tiled.range)
}

makeTiledRange <- function(gr, binwidth){
  # Bin Range, currently only supports 1 range
  tiled.range <- tile(gr, width=binwidth)[[1]]
  
  # Remove strand information
  strand(tiled.range) <- "*"
  
  return(tiled.range)
}

binCoverageInRange <- function(cvg.L, gr, binwidth=1000, val="score", 
                               scaling=log){
  # Bin Range, currently only supports 1 range
  tiled.range <- tile(gr, width=binwidth)[[1]]
  
  # Remove strand information
  strand(tiled.range) <- "*"
  
  for(g in names(cvg.L)){
    # Find overlap between the coverage and the tiled range
    overlap <- findOverlaps(cvg.L[[g]], tiled.range)
    
    # Split scores from coverage range into their appropriate bin and sum
    bin.split <- splitAsList(mcols(cvg.L[[g]])[[val]][queryHits(overlap)],
                             factor(subjectHits(overlap)))
    
    # Need a better way to scale the data
    bin.score <- scaling(sum(bin.split))
    elementMetadata(tiled.range)[[g]] <- bin.score
  }
  return(tiled.range)
}

#' makeDataTracks
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
makeDataTracks <- function(cvg.L, gr, genome, chr, bg.title, score.max, colors = NULL,
                           type=NULL, binwidth=1000, val="score", scaling=log,
                           dtrack.kwargs=list()){
  # Compile sample coverages as datatracks
  dtrack <- list()
  if(type == 'heatmap'){
    hm.gradient <- colorRampPalette(brewer.pal(9, "BuGn"))(100)
    binned.cvg <- binCoverageInRange(cvg.L, gr, binwidth, val, scaling)
    heatmap.params <- list(range=binned.cvg, type='heatmap', chromosome=chr,
                           background.title='gray50', showSampleNames = TRUE,
                           cex.sampleNames = 0.6, cex.axis=0.6, gradient = hm.gradient)
    heatmap.params <- modifyList(heatmap.params, dtrack.kwargs)
    dtrack[[1]] <- do.call(Gviz::DataTrack, heatmap.params)
  } else if(type == 'hist'){
    for (g in names(cvg.L)) {
      coverage.params <- list(start=start(cvg.L[[g]]), end=end(cvg.L[[g]]),
                              data=score(cvg.L[[g]]), chromosome=chr, type='hist',
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