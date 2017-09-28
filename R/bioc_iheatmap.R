
# makes x based on colnames of mat if available
# if not available, just uses 1 to number of columns
se_default_x <- function(se, assay){
  mat <- assays(se)[[assay]]
  if (is.null(colnames(mat))){
    return(seq_len(ncol(mat)))
  } else{
    colnames(mat)
  }
}
# makes y based on rownames of mat if available
# if not available, just uses 1 to number of rows
se_default_y <- function(se, assay){
  mat <- assays(se)[[assay]]
  if (is.null(rownames(mat))){
    return(seq_len(nrow(mat)))
  } else{
    rownames(mat)
  }
}

#' Summarized Experiment iheatmap and iheatmap methods
#' 
#' @param p IHeatmap object
#' @param data SummarizedExperiment
#' @param assay name of assay in data to plot
#' @param x name of rows of data, by default rownames
#' @param y name of columns of data, by default colnames
#' @param cluster_rows how to cluster rows, default is hclust
#' @param cluster_cols how to cluster columns, default is hclust
#' @param scale scale rows, columns, both, or none
#' @param ... additional arguments to \code{\link[iheatmapr]{iheatmap}}
#' @rdname iheatmap-SummarizedExperiment
#' @name iheatmap-SummarizedExperiment
#' @aliases iheatmap,SummarizedExperiment-method
#' add_iheatmap,IHeatmapHorizontal,SummarizedExperiment-method
#' add_iheatmap,IHeatmapVertical,SummarizedExperiment-method
#' @export
setMethod("iheatmap", c(data = "SummarizedExperiment"),
          function(data, 
                   assay = assayNames(data)[[1]],
                   x = se_default_x(data, assay),
                   y = se_default_y(data, assay),                   
                   cluster_rows = c("hclust","kmeans","none"),
                   cluster_cols = c("hclust","kmeans","none"),
                   scale = "rows",
                   ...){

            stopifnot(length(assay) == 1)
            stopifnot(assay %in% assayNames(data))
            mat <- assays(data)[[assay]]
            iheatmap(mat, x, y, cluster_rows = match.arg(cluster_rows),
                            cluster_cols = match.arg(cluster_cols),
                            scale = "rows",
                     ...)  
         
          })

#' @rdname iheatmap-SummarizedExperiment
#' @export
setMethod("add_iheatmap", c(p = "IheatmapHorizontal", 
                            data = "SummarizedExperiment"),
          function(p, data, 
                   assay = assayNames(data)[[1]],
                   x = se_default_x(data),
                   cluster_cols = c("hclust","kmeans","none"),
                   scale = "rows",
                   ...){
            
            stopifnot(length(assay) == 1)
            stopifnot(assay %in% assayNames(data))
            mat <- assays(data)[[assay]]
            add_iheatmap(p, mat, x, 
                     cluster_cols = match.arg(cluster_cols),
                     scale = scale, ...)  
            
          })

#' @rdname iheatmap-SummarizedExperiment
#' @export
setMethod("add_iheatmap", c(p = "IheatmapVertical", 
                            data = "SummarizedExperiment"),
          function(p,
                   data, 
                   assay = assayNames(data)[[1]],
                   y = se_default_y(data),
                   cluster_rows = c("hclust","kmeans","none"),
                   scale = "rows",
                   ...){
            
            stopifnot(length(assay) == 1)
            stopifnot(assay %in% assayNames(data))
            mat <- assays(data)[[assay]]
            add_iheatmap(p, mat, y, cluster_rows = match.arg(cluster_rows),
                      scale = scale,
                     ...)  
            
          })
