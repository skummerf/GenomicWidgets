##' Expand sample table with paths to useful files
##'
##' Expand sample table with paths to useful files
##' @title Expand sample table with paths to useful files
##' @param dat.samples Input sample table with minimal columns Sample.Name and Dir.ngs_pipeline
##' @return data frame including additional columns
##' @author Suchit Jhunjhunwala
##' @export
expandSampleTable <- function(dat.samples, dir.bed) {
  requiredColumns <- c("Sample.Name", "Dir.ngs_pipeline")
  if (! all(requiredColumns %in% colnames(dat.samples))) {
    logerror(paste("Required columns missing from sample table:",
                   setdiff(requiredColumns, colnames(dat.samples))
    )
    )
    return(NULL)
  }
  
  if (! missing(dir.bed) ) {
    dat.samples$File.reads.bed <- paste0(dir.bed,'/',dat.samples$Sample.Name,'.bed')
  }
  
  dat.samples <- within(dat.samples, {
    if (! "File.cvg" %in% colnames(dat.samples)) File.cvg <- rep(NA, nrow(dat.samples))
    File.cvg <- ifelse(is.na(File.cvg),
                       paste0(Dir.ngs_pipeline, "/results/", Sample.Name, ".coverage.RData"),
                       File.cvg
    )
    
    if (! "File.bam" %in% colnames(dat.samples)) File.bam <- rep(NA, nrow(dat.samples))
    File.bam <- ifelse(is.na(File.bam),
                       paste0(Dir.ngs_pipeline, "/bams/", Sample.Name, ".analyzed.bam"),
                       File.bam
    )
    if (! "File.bw" %in% colnames(dat.samples)) File.bw <- rep(NA, nrow(dat.samples))
    File.bw <- ifelse(is.na(File.bw),
                      paste0(Dir.ngs_pipeline, "/results/", Sample.Name, ".coverage.bw"),
                      File.bw
    )
    
    if (! "File.normalizedBw" %in% colnames(dat.samples)) File.normalizedBw <- rep(NA, nrow(dat.samples))
    File.normalizedBw <- ifelse(is.na(File.normalizedBw),
                                paste0(Dir.ngs_pipeline, "/results/", Sample.Name, ".coverage.normalized.bw"),
                                File.normalizedBw
    )
    
    
  })
  
  return(dat.samples)
}


##' expand the paired sample table with default directory structure
##'
##' expand the paired sample table with default directory structure, adding directories for
##' macs output, peak summary, and peak motif. Directories predifined by the user will
##' be preserved
##' @param dat.pair input paired sample table
##' @param order sample.first means the top directory will be named by the sample
##' pair, subdir names are macs_output etc.
##' @return data frame of expanded sample table
##' @export
##' @author Jinfeng Liu
expandPairedTable <- function(dat.pair, order='sample.first') {
  dirNames <- list(macs='macs_output',
                   peakSummary='peak_summary',
                   peakMotif='peak_motif')
  
  for (d in names(dirNames)) {
    col.name <- paste0('Dir.',d)
    if ( ! col.name %in% colnames(dat.pair) ) {
      dat.pair[,col.name] <- NA
    }
    
    for (i in 1:nrow(dat.pair) ) {
      # skip rows with user-defined dir or Dir.chipseq with NA
      if ( ! is.na(dat.pair[i,col.name]) | is.na(dat.pair[i,'Dir.chipseq']) |
           is.na(dat.pair[i,'Prefix.pairname']))
        next
      
      if ( order == 'sample.first' ) {
        dat.pair[i,col.name] <-
          paste0(dat.pair[i,"Dir.chipseq"],'/',
                 dat.pair[i,'Prefix.pairname'],'/',dirNames[[d]],'/')
      } else {
        dat.pair[i,col.name] <-
          paste0(dat.pair[i,"Dir.chipseq"],'/',
                 dirNames[[d]],'/',dat.pair[i,'Prefix.pairname'],'/')
      }
    }
  }
  
  dat.pair <- within(dat.pair, {
    if(! "File.macs" %in% colnames(dat.pair)) {
      File.macs <- rep(NA, nrow(dat.pair))
    }
    File.macs <-
      ifelse(is.na(File.macs) & ! is.na(Prefix.pairname),
             file.path(Dir.macs, paste0(Prefix.pairname, "_peaks.xls")),
             File.macs
      )
    
    if(! "File.peaks" %in% colnames(dat.pair)) {
      File.peaks <- rep(NA, nrow(dat.pair))
    }
    File.peaks <-
      ifelse(is.na(File.peaks) & ! is.na(Prefix.pairname),
             file.path(Dir.peakSummary, "annotated_peaks.txt"),
             File.peaks
      )
    
  })
  
  
  return(dat.pair)
}

##' rename MACS output column names
##'
##' rename MACS output column names
##' @param x a vector of MACS column names
##' @return a vector of fixed column names
##' @export 
##' @author Jinfeng Liu
renameMacsColname <- function(x) {
  x[x == "-10*log10(pvalue)"] <- 'neg10.log10P'
  x[x == "FDR(%)"] <- "FDR"
  x[x == 'tags' ] <- 'pileup'
  x[x == 'fold_enrichment'] <- 'fold.enrichment'
  x[x == '-log10(qvalue)'] <- 'neg.log10Q'
  x[x == '-log10(pvalue)'] <- 'neg.log10P'
  return(x)
}


##' Read peaks from MACS output
##'
##' 
##' Read peaks from MACS output, both v1.4 and v2. Column names in v1 and v2 will be harmonized, use the peak start/end position as the main positions by default
##' @param file filename
##' @param fc Fold-enrichment threshold 
##' @param minTags Minimum number of tags in the peak
##' @param neglog10pThreshold Threshold value for the "-log10(p.value)" reported in MACS output
##' @param fdr FDR cut-off
##' @param use.summit use summit position instead of peak start/end position as main positions
##' @param broad whether the input is macs2's broadpeak file
##' @param as.Granges Logical. If the peaks should be returned as a GRanges object.
##' @return data frame of macs peaks, with cleaned up column names and FDR values converted from percentage to relative value (between 0 and 1)
##' @author Suchit Jhunjhunwala, Jinfeng Liu
##' @export
##' @import GenomicRanges
readMacsPeaks <- function(file, fc, minTags, neglog10pThreshold, fdr,
                          use.summit=FALSE,as.Granges=TRUE,broad=FALSE) {
  
  all.chrom <- c(1:22,'X','Y')
  if ( ! broad ) {
    df <- tryCatch({
      read.table(file,as.is=TRUE,header=TRUE,sep="\t",check.names=FALSE)
    }, error = function(e) {
      print(sprintf('Error reading file %s returning NULL', file))
      return(NULL)
    })
  } else {
    df <- tryCatch({
      read.table(file,as.is=TRUE,header=FALSE,sep="\t")
    }, error = function(e) { 
      print(sprintf('Error reading file %s returning NULL', file))
      return(NULL)
    })
    colnames(df) <- c("chr","start","end","name","value","notused","fold.enrichment",
                      "neg.log10P","neg.log10Q")
    df <- df[,setdiff(colnames(df),c('value','notused'))]
  }
  
  if(is.null(df))
    return(df)
  
  df$chr <- sub("chr","",df$chr)
  df <- df[df$chr %in% all.chrom, ]
  
  colnames(df) <- renameMacsColname(colnames(df))
  
  if ( 'FDR' %in% colnames(df) ) {
    df$FDR <- df$FDR / 100
    df$neg.log10Q <- -log10(df$FDR)
  } else {
    df$FDR <- 10^(-1 * df$neg.log10Q)
  }
  
  
  if ( ! 'neg.log10P' %in% colnames(df)) {
    df$neg.log10P <- df$neg10.log10P / 10
  }
  if ( ! broad ) {
    if ( 'summit' %in% colnames(df) ) {
      df$summit <- df$start + df$summit
    } else {
      df$summit <- df$abs_summit
    }
  }
  
  if ( ! 'pileup' %in% colnames(df) && 'tags' %in% colnames(df) ) {
    df$pileup <- df$tags
  }
  
  
  if(! missing(fdr)) {
    df <- df[df$FDR <= fdr, ]
  }
  if(! missing(fc)) {
    df <- df[df$fold.enrichment >= fc, ]
  }
  if(! missing(minTags) && 'pileup' %in% colnames(df) ) {
    df <- df[ df$pileup >= minTags, ]
  }
  if(! missing(neglog10pThreshold)) {
    df <- df[df$neg.log10P >= neglog10pThreshold, ]
  }
  
  if (! 'name' %in% colnames(df)) { 
    df$name <- paste0("peak", seq_len(nrow(df)))
  }
  
  rownames(df) <- df$name
  df$peak.start <- df$start
  df$peak.end <- df$end
  
  if (! 'length' %in% colnames(df)) df$length <- df$end - df$start + 1
  
  if ( ! broad && use.summit && 'summit' %in% colnames(df) ) {
    df$start <- df$summit
    df$end <- df$summit
  }
  
  columns.out <- c("chr","start","end","length","summit","peak.start","peak.end",
                   "pileup","neg.log10P","neg.log10Q","FDR","fold.enrichment","name")
  
  if ( broad ) {
    columns.out <- setdiff(columns.out, c("peak.start","peak.end","pileup"))
  }
  
  if ( ! 'summit' %in% colnames(df)) {
    columns.out <- setdiff(columns.out, 'summit')
  }
  
  df <- df[,columns.out]
  
  
  if(as.Granges==FALSE) {
    return(df)
  } else {
    gr <- with(df, GRanges(
      seqnames=chr,
      ranges=IRanges(start, end)
    ))
    mcols(gr) <- df[, -1:-3]
    return(gr)
  }
}


##Function from gChipseq
ngsPipelineResult <- function (dirs, file = ".summary_alignment.tab$", sep = "\t", 
                               column = "value", brief.name = TRUE, ncores = 1) {
  ncores <- ifelse(is.null(ncores), detectCores(), min(ncores, 
                                                       detectCores()))
  value.list <- list()
  dirs <- unique(dirs)
  value.list <- mclapply(setNames(dirs, dirs), function(dir) {
    dir.result <- file.path(dir, "results")
    file.in <- list.files(dir.result, pattern = file)
    if (length(file.in) != 1) {
      logerror("number of files not equal to 1")
      return(NULL)
    }
    file.in <- file.path(dir.result, file.in)
    dat.in <- read.table(file.in, as.is = TRUE, header = TRUE, 
                         row.names = 1, sep = "\t")
    if (nrow(dat.in) < 1) 
      next
    if (!column %in% colnames(dat.in)) {
      logerror(paste(column, " not in the column of", file, 
                     "\n"))
      return(NULL)
    }
    dat.out <- dat.in[, column]
    names(dat.out) <- rownames(dat.in)
    return(dat.out)
  }, mc.cores = ncores)
  if (length(value.list) < 1) {
    logerror("no data read in")
    return(NULL)
  }
  nrow.data <- unique(sapply(value.list, length))
  rownames.data <- names(value.list[[1]])
  if (length(nrow.data) != 1) {
    logerror("inconsistent row numbers")
    return(NULL)
  }
  dat.out <- matrix(NA, nrow = nrow.data, ncol = length(value.list), 
                    dimnames = list(rownames.data, names(value.list)))
  for (i in 1:length(value.list)) {
    dat.out[, i] <- value.list[[i]][rownames.data]
  }
  if (brief.name) {
    samples <- colnames(dat.out)
    samples <- sub("\\/$", "", samples)
    samples <- sub(".*\\/", "", samples)
    colnames(dat.out) <- samples
  }
  if (ncol(dat.out) != length(dirs)) {
    logwarn(paste(length(dirs), "input dirs,", ncol(dat.out), 
                  " output columns\n"))
  }
  return(dat.out)
}

##' get tss from txdb
##'
##' get tss information from txdb
##' @title get tss information from txdb
##' @param txdb input txdb object
##' @param regularize Logical, indicating whether only the regular chromosomes should be output
##' @return tss granges
##' @author Jinfeng Liu, Suchit Jhunjhunwala
##' @import GenomicFeatures
##' @export
tssFromTxdb <- function(txdb, regularize=TRUE) {
  tx <- as.data.frame(transcripts(txdb))
  tx$seqnames <- sub("chr","",tx$seqnames)
  if(regularize) {
    regular_chr <- c(1:26,"X","Y")
    tx <- tx[tx$seqnames %in% regular_chr, ]
  }
  tss <- GRanges(seqnames=tx$seqnames,
                 ranges=IRanges(
                   start=ifelse(tx$strand=="+", tx$start, tx$end),
                   width=1),
                 names = tx$tx_name, strand =tx$strand)
  return(tss)
}