
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
plotGeneCoverage <- function(grl,symbol,org,genome,extend=0,
                             col='blue',bg.title='black',
                             sync=FALSE,ymax,scale.group=1,
                             dat.exon, ...) {
  
  # Check that the genome is available
  if(!speciesGenomeMatch(org, genome))
    return(NULL)
  
  # Select the correct genome
  mart <- ifelse (genome %in% c("GRCh38","GRCm38"), "igis", "igis2.3")
  
  # Import exon data if missing
  if ( missing(dat.exon) ) {
    dat.exon <- igisExonlist(symbol,org=org, mart=mart, ...)
  }
  
  # Exit if no exon data was provided or imported
  if ( nrow(dat.exon) == 0 ) {
    warning(paste("no exon data for",symbol,"\n"))
    return(NULL)
  }
  
  # Make the GRange object for the gene
  range.gene <- GRanges(seqnames=dat.exon[1,"chr"],
                        ranges=IRanges(start=min(dat.exon$start),
                                       end=max(dat.exon$end)),
                        strand=dat.exon[1,"strand"])
  
  # range.gene <- grAddChr(range.gene)
  
  # Ensure that the extension is a length 2 vector
  if (length(extend) == 1) {
    extend <- rep(extend,2)
  } else if ( length(extend) > 2) {
    warning("extend should be a vector of length 1 or 2")
    return(NULL)
  }
  
  # Add the extension to the gene range for each strand
  if (as.character(strand(range.gene)) == '+') {
    start(range.gene) <- start(range.gene) - extend[1]
    end(range.gene) <- end(range.gene) + extend[2]
  } else if ( as.character(strand(range.gene)) == '-' ) {
    end(range.gene) <- end(range.gene) + extend[1]
    start(range.gene) <- start(range.gene) - extend[2]
  }
  
  # Get the chromosome
  chr <- as.character(seqnames(range.gene))
  
  # Make colors for the plots
  colors <- rep(col,length.out=length(grl))
  names(colors) <- names(grl)
  
  # Get the gene coverage for each sample supplied
  gene.cov <- list()
  for (g in names(grl) ) {
    # grl[[g]] <- grAddChr(grl[[g]])
    gene.overlap <- findOverlaps(query=range.gene, grl[[g]])
    gene.cov[[g]] <- grl[[g]][as.data.frame(gene.overlap)[,2],]
  }
  
  # Scale data range
  if ( missing(ymax) ) {
    score.max <- sapply(gene.cov, function(x) { max(score(x)) } )
    if ( sync ) {
      scale.group <- rep(scale.group, length.out=length(gene.cov))
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
  
  # Compile sample coverages as datatracks
  dtrack <- list()
  for (g in names(grl)) {
    dtrack[[g]] <- DataTrack(start=start(gene.cov[[g]]),end=end(gene.cov[[g]]),
                             data=score(gene.cov[[g]]), chromosome=chr, type='hist',
                             ylim=c(0,score.max[g]),background.title=bg.title[g],
                             genome=genome, name=g,col.histogram=colors[g],
                             fill.histogram=colors[g])
  }
  
  # Add genome tracks
  gtrack <- GenomeAxisTrack()
  grtrack <- GeneRegionTrack(dat.exon, genome=genome,chromosome=chr,
                             name=symbol,background.title='orange')
  # Add tracks
  tracklist <- list()
  tracklist <- c(gtrack,grtrack)
  for (g in names(grl)) {
    tracklist <- c(tracklist,dtrack[[g]])
  }
  
  # Plot tracks
  plotTracks(tracklist,main=symbol,
             from=start(range.gene),to=end(range.gene))
  
}

# Make track

# 1. Genome track
# 2. Annotation Track
# 3. Chip Coverage
# 4. RNA Coverage

