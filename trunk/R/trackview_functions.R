
##' plot genomic coverage along the gene
##' Adapted from gChipseq function of the same name. Credit to Jinfeng Liu
##' plot genomic coverage along the gene
##' @param grl a list of GRangesList, typically the results of importing wig/bw/... files for each sample
##' @param symbol str; symbol to plot
##' @param org str; organism, human|mouse
##' @param genome str; genome build, e.g. hg19, GRCh38, etc.
##' @param extend extend the plotting region
##' @param col colors of the coverage plot
##' @param bg.title background color for the title panel
##' @param sync logical for whether to sync the ymax for all data tracks
##' @param ymax vector ymax for the data tracks
##' @param scale.group groups to use the same y scale
##' @param dat.exon a data frame representing exon information, if NULL get it from IGIS
##' @param ... 
##' @return nothing
##' @import Gviz gChipseq
##' @import IRanges
##' @export
##' @author Justin Finkle
plotGeneCoverage <- function(grl, dat.exon, target.range, symbol, genome, 
                             col='blue', ymax, bg.title='black', sync=FALSE, 
                             scale.group=1, hm.thresh=4) {
  
  # Get the chromosome
  chr <- as.character(seqnames(target.range))
  
  # Make colors for the plots
  colors <- rep(col,length.out=length(grl))
  names(colors) <- names(grl)
  
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
  
  dtrack <- makeDataTracks(grl, target.range, genome, chr, bg.title, colors, 
                           score.max = score.max, hm.thresh=hm.thresh)
  

  
  # Add genome tracks
  grtrack.color <- 'lightslateblue'
  gtrack <- GenomeAxisTrack()
  target.range <- grAddChr(target.range)
  options(ucscChromosomeNames = FALSE)
  if (length(dtrack)==1){
    ann.size <- 0.25
  } else {
    ann.size <- NULL
  }
  grtrack <- GeneRegionTrack(dat.exon, genome=genome, chromosome=chr, name=symbol,
                             background.title=grtrack.color, fill=grtrack.color,
                             transcriptAnnotation="symbol", size=ann.size)
  # Add tracks
  tracklist <- list()
  tracklist <- c(gtrack,grtrack, dtrack)
  # if (length(dtrack) >1){
  #   for (g in names(grl)) {
  #     tracklist <- c(tracklist,dtrack[[g]])
  #   }
  # }
  # Plot tracks
  ptype <- ifelse(length(grl) > hm.thresh, 'heatmap', 'histogram')
  plotTracks(tracklist,main=symbol, from=start(target.range),
             to=end(target.range), showSampleNames = TRUE, 
              cex.sampleNames = 0.6)
  
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
getCoverageInRange <- function(bwList, target.range, names){
  # Import only the range that matches target.ange
  cvg <- lapply(bwList, function(x) import.bw(x, which=target.range))
  cvg <- GRangesList(cvg)
  
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


binCoverageInRange <- function(cvg.L, gr, binwidth=1000, val="score", 
                               scaling=log){
  # Bin Range, currently only supports 1 range
  tiled.range <- tile(gr, width=binwidth)[[1]]
  
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

makeDataTracks <- function(cvg.L, gr, genome, chr, bg.title, colors,score.max,
                           hm.thresh=4, binwidth=1000, val="score",
                           scaling=log, ...){
  # Compile sample coverages as datatracks
  dtrack <- list()
  if(length(cvg.L) > hm.thresh){
    binned.cvg <- binCoverageInRange(cvg.L, gr, binwidth, val, scaling)
    dtrack[[1]] <- DataTrack(binned.cvg, type='heatmap', chromosome=chr)
  } else {
    for (g in names(cvg.L)) {
      dtrack[[g]] <- DataTrack(start=start(cvg.L[[g]]),end=end(cvg.L[[g]]),
                               data=score(cvg.L[[g]]), chromosome=chr, type='hist',
                               ylim=c(0,score.max[g]),background.title=bg.title[g],
                               genome=genome, name=g, col.histogram=colors[g],
                               fill.histogram=colors[g]) 
    }
    names(dtrack) = names(cvg.L)
  }
  return(dtrack)
}
