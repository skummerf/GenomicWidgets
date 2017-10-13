setClass("TranscriptParts",
         slots = c(intron = "GRangesList",
          utr5 = "GRangesList",
          utr3 = "GRangesList",
          cds = "GRangesList",
          exon = "GRangesList",
          transcripts = "GRanges",
          seqlevelsStyle =  "character")
)


setClass("ViewRange",
         slots = c(range = "GRanges",
                   relative = "logical",
                   reference = "numeric"))

setClass("LocusPlot",
         slots = c(layout = "list",
                   trackname = "character"
         ),
         contains = "VIRTUAL")


setClass("AnnotationPlot",
         slots = c(transcripts = "GRanges"),
         contains = "LocusPlot")

setClass("SignalPlot",
         slots = c(signal = "SimpleList",
                   x = "numeric",
                   mode = "character",
                   fill = "character",
                   color = "character",
                   showlegend = "logical"
         ),
         contains = "LocusPlot")

setClass("LocusSummary",
         slots = c(data = "list",
                   layout = "list"))

#' @importClassesFrom S4Vectors SimpleList
setClass("LocusView",
         slots = c(view = "ViewRange",
                   layout = "list",
                   share_y = "logical",
                   heights = "numeric"),
         prototype = list(elementType = "LocusPlot"),
         contains = c("SimpleList"))


setClass("LocusViewList",
         slots = c(share_y = "logical",
                   xtitle = "character"),
         prototype = list(elementType = "LocusView"),
         contains = c("SimpleList"))

setClass("LocusSummaryList",
         prototype = list(elementType = "LocusSummary"),
         contains = c("SimpleList"))

setClassUnion("TranscriptPartsOrNull", c("TranscriptParts", "NULL"))



setClass("SummaryParameters",
         slots = list(
           data = "SummarizedExperiment",
           assay_name = "character",
           groups = "factor",
           showlegend = "logical",
           colors = "character",
           boxpoints = "character",
           pointpos = "numeric",
           ytitle = "character",
           width = "numeric",
           ranges = "GenomicRanges"
         ))

setClassUnion("SummaryParametersOrNull", c("SummaryParameters", "NULL"))

setClass("TrackParameters",
         slots = list(
           data = "character",
           annotation = "TranscriptPartsOrNull",
           track_names = "character",
           groups = "factor",
           share_y = "logical",
           showlegend = "logical",
           colors = "character",
           fill = "character",
           mode = "character",
           annotation_position = "character",
           annotation_size = "numeric",
           summary = "SummaryParametersOrNull",
           layout = "list"
         ))






setClass("GenomeTrackWidget",
         slots = c(tracks = "LocusViewList",
                   summaries = "LocusSummaryList",
                   summary_width = "numeric",
                   layout = "list"))




