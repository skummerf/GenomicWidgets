# SummarizedExperiment Method

setMethod("make_locus_summary", c("SummarizedExperiment"),
          function(object,  assay = assayNames(object)[1],
                   ..., 
                   groups = NULL,
                   type = c("auto","box","bar","dot","density"),
                   showlegend = TRUE, 
                   colors = NULL,
                   signal = "Expression"){
            
            # Make LocusSummary
            new("SignalPlot",
                signal = sm,
                color = colors,
                mode = mode,
                fill = fill,
                trackname = name)
          })



# eset / ExpressionSet Method
