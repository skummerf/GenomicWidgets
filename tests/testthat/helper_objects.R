#data("rpkm_chr21", package = "GenomicWidgets")

genomation_dir <- system.file("extdata", package = "genomationData")

samp.file <- file.path(genomation_dir,'SamplesInfo.txt')
samp.info <- read.table(samp.file, header=TRUE, sep='\t',
                        stringsAsFactors = FALSE)
samp.info$fileName <- file.path(genomation_dir, samp.info$fileName)

ctcf.peaks <- genomation::readBroadPeak(
  system.file("extdata",
              "wgEncodeBroadHistoneH1hescCtcfStdPk.broadPeak.gz",
              package = "genomationData"))
ctcf.peaks <- ctcf.peaks[GenomicRanges::seqnames(ctcf.peaks) == "chr21"]
ctcf.peaks <- ctcf.peaks[order(-ctcf.peaks$signalValue)]
ctcf.peaks <- GenomicRanges::resize(ctcf.peaks, width = 501, fix = "center")

#tss <- GenomicRanges::promoters(SummarizedExperiment::rowRanges(rpkm_chr21), 
#                                 up = 1, down = 1)