setClass("TranscriptParts",
         slots = c(intron = "GRangesList",
          utr5 = "GRangesList",
          utr3 = "GRangesList",
          cds = "GRangesList",
          exon = "GRangesList",
          transcripts = "GRanges",
          seqlevelsStyle =  "character")
)

#' @importClassesFrom GenomicRanges GRanges
setClass("ViewRange",
         contains = "GRanges")

setValidity("ViewRange", 
            function(object){
              "name" %in% colnames(mcols(object))
            })

setClass("RelativeViewRange",
         contains = "ViewRange")

setValidity("RelativeViewRange", 
             function(object){
               "relative" %in% colnames(mcols(object))
             })


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
                   mode = "character",
                   fill = "character",
                   color = "character",
                   showlegend = "logical"
         ),
         contains = "LocusPlot")

setClass("LocusSummary",
         slots = c(data = "list",
                   layout = "list"))

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

setClass("GenomeTrackWidget",
         slots = c(tracks = "LocusViewList",
                   summaries = "LocusSummaryList",
                   summary_width = "numeric",
                   layout = list()))




