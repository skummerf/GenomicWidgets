% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/track_interface.R
\docType{methods}
\name{set_summary_parameters}
\alias{set_summary_parameters}
\alias{set_summary_parameters,SummarizedExperiment-method}
\alias{set_summary_parameters,RangedSummarizedExperiment-method}
\alias{set_summary_parameters,RangedSummarizedExperiment-method}
\title{set_summary_parameters}
\usage{
\S4method{set_summary_parameters}{SummarizedExperiment}(object, ranges,
  assay_name = assayNames(object)[1], groups = colnames(object),
  colors = "blue", showlegend = length(colors) > 1, boxpoints = c("all",
  "Outliers", "false"), pointpos = 0, ytitle = "Expression", width = 0.3)

\S4method{set_summary_parameters}{RangedSummarizedExperiment}(object,
  ranges = rowRanges(object), assay_name = assayNames(object)[1],
  groups = NULL, colors = "blue", showlegend = length(colors) > 1,
  boxpoints = c("all", "Outliers", "false"), pointpos = 0,
  ytitle = "Expression", width = 0.3)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{ranges}{ranges corresponding to rows of object}

\item{assay_name}{name of assay to use}

\item{groups}{either vector of group assignments of name of column in object 
colData that corresponds to vector of group assignments}

\item{colors}{colors to use}

\item{showlegend}{show the legend?}

\item{boxpoints}{plot individual points?}

\item{pointpos}{relative position of points to boxes}

\item{ytitle}{name for yaxis}

\item{width}{relative width of summary plots when plotting tracks}

\item{...}{additional arguments}
}
\value{
object storing summary parameters, for use in
 \code{\link{set_track_parameters}}
}
\description{
Setup parameters for plotting summaries along genome tracks.
Result from this function can be passed to \code{\link{set_track_parameters}}
}
\examples{


library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

## we'll read in some RNA counts
data(rpkm_chr21)

## From the ranges of the rpkm object, we'll pull out the tss
chr21_promoters <- promoters(SummarizedExperiment::rowRanges(rpkm_chr21),
                 up = 1000, down = 1000)

## set summary parameters

summary_params <- set_summary_parameters(rpkm_chr21,
  groups = "GROUP", ranges = chr21_promoters) 
  
## We'll also read in some track data to plot
genomation_dir <- system.file("extdata", package = "genomationData")
samp.file <- file.path(genomation_dir,'SamplesInfo.txt')
samp.info <- read.table(samp.file, header=TRUE, sep="\\t", 
                        stringsAsFactors = FALSE)
samp.info$fileName <- file.path(genomation_dir, samp.info$fileName)

## Make track plotter using summary parametrs

track_params <- set_track_parameters(samp.info$fileName[1:3], 
  annotation = TxDb.Hsapiens.UCSC.hg19.knownGene, 
  track_names = samp.info$sampleName[1:3], 
  share_y = TRUE,
  summary = summary_params)
  
if (interactive()){
  plot_tracks(rownames(rpkm_chr21)[1:3], track_params)
  plot(chr21_promoters[1:3], track_params)
}   

}
\author{
Alicia Schep and Justin Finkle
}
