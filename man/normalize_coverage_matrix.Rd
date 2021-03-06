% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coverage_matrix.R
\docType{methods}
\name{normalize_coverage_matrix}
\alias{normalize_coverage_matrix}
\alias{normalize_coverage_matrix,list-method}
\alias{normalize_coverage_matrix,matrix-method}
\alias{normalize_coverage_matrix,SummarizedExperiment-method}
\alias{normalize_coverage_matrix,matrix-method}
\alias{normalize_coverage_matrix,SummarizedExperiment-method}
\title{normalize_coverage_matrix}
\usage{
\S4method{normalize_coverage_matrix}{list}(mats, method = c("localRms",
  "localMean", "localNonZeroMean", "PercentileMax", "scalar", "none"),
  pct = 0.95, scalar = NULL, digits = 3)

\S4method{normalize_coverage_matrix}{matrix}(mats, method = c("localRms",
  "localMean", "localNonZeroMean", "PercentileMax", "scalar", "none"),
  pct = 0.95, scalar = NULL, digits = 3)

\S4method{normalize_coverage_matrix}{SummarizedExperiment}(mats, ...)
}
\arguments{
\item{mats}{matrix, list of matrix, or SummarizedExperiment}

\item{method}{normalization method option, see Details}

\item{pct}{Percentile, only used if PercentileMax is method}

\item{scalar}{vector of scalars used for normalizing each mat, only 
used if scalar is method}

\item{digits}{number of significant digits of result to keep.}

\item{...}{additional arguments to normalize_coverage_matrix}
}
\value{
Should return data in the same format as input, but now with values
normalized according to the method chosen.
}
\description{
Normalizes coverage matrices using one of several methods.
}
\details{
Normalization choices are "localRms", "localMean", 
"localNonZeroMean", "PercentileMax", "scalar", and "none".  localRMS will 
divide each row by the root mean squared values of that row.  localMean will
divide each row by the mean of that row.  localNonZeroMean will divide each 
row by nonzero values in that row.  PercentileMax will divide values based on 
percentile (given by pct argument) of the entire matrix.  scalar will divide
entire matrix by a scalar, given by scalar argument.  This scalar could for 
example be a measure of the sequencing depth.
}
\examples{
## First we'll make some coverage matrices

library(GenomicRanges)
# First read in some sample data
genomation_dir <- system.file("extdata", package = "genomationData")

samp.file <- file.path(genomation_dir,'SamplesInfo.txt')
samp.info <- read.table(samp.file, header=TRUE, sep='\\t', 
                       stringsAsFactors = FALSE)
samp.info$fileName <- file.path(genomation_dir, samp.info$fileName)

ctcf.peaks = genomation::readBroadPeak(system.file("extdata",
               "wgEncodeBroadHistoneH1hescCtcfStdPk.broadPeak.gz",
               package = "genomationData"))
ctcf.peaks = ctcf.peaks[seqnames(ctcf.peaks) == "chr21"]
ctcf.peaks = ctcf.peaks[order(-ctcf.peaks$signalValue)]
ctcf.peaks = resize(ctcf.peaks, width = 1000, fix = "center")

# Make the coverage matrices
mats <- make_coverage_matrix(samp.info$fileName[1:3], ctcf.peaks, 
                     up = 500, down = 500, binsize = 25)
                     
# Now normalize:
norm_mats <- normalize_coverage_matrix(mats)                    
                     
}
\author{
Alicia Schep
}
