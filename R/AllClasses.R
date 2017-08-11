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
         contains = "GRanges")

setValidity("ViewRange", 
            function(object){
              "name" %in% colnames(mcols(object))
            })

setClass("RelativeViewRange",
         contains = "ViewRange")

setValidity("RelativeViewRange", 
            function(object){
              ("relative" %in% colnames(mcols(object)) & "reference" %in% 
                 colnames(mcols(object)))
            })


setClass("LocusPlot",
         slots = c(layout = "list",
                   trackname = "character"
         ),
         contains = "VIRTUAL")


setClass("AnnotationPlot",
         slots = c(transcripts = "GRanges"),
         contains = "LocusPlot")

#' @importClassesFrom genomation ScoreMatrixList
setClass("SignalPlot",
         slots = c(signal = "ScoreMatrixList",
                   mode = "character",
                   fill = "character",
                   color = "character",
                   showlegend = "logical"
         ),
         contains = "LocusPlot")

setClass("LocusSummaryPlot",
         slots = c(data = "list",
                   layout = "list"))

setClass("LocusView",
         slots = c(view = "ViewRange",
                   layout = "list",
                   share_y = "logical",
                   summary = "LocusSummaryPlot",
                   summary_width = "numeric",
                   heights = "numeric"),
         prototype = list(elementType = "LocusPlot"),
         contains = c("SimpleList"))



setClass("MultiLocusView",
         slots = c(layout = "list",
                   share_y = "logical"),
         prototype = list(elementType = "LocusView"),
         contains = c("SimpleList"))




