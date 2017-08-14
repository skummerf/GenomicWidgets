# SummarizedExperiment Method

setMethod("make_locus_summary", c("SummarizedExperiment"),
          function(object,  
                   row_name,
                   assay_name = assayNames(object)[1],
                   ..., 
                   groups = NULL,
                   showlegend = !is.null(colors), 
                   colors = NULL,
                   boxpoints = c("all","Outliers","false"),
                   pointpos = 0,
                   ytitle = "Expression"){
            
            boxpoints = match.arg(boxpoints)
            
            # Check groups
            if (is.null(groups)){
              groups <- ""
            } else if (length(groups) == 1){
              if (groups %in% colnames(colData(object))) groups <- colData(object)[,groups]
            } else{
              stopifnot(length(groups) == ncol(object))
            }
            
            
            if (is.null(colors) || length(colors) == 1){
              data <- list(list(y = assay(object,assay_name)[row_name,],
                         x = groups,
                         type = "box",
                         text = text,
                         name = ytitle,
                         pointpos = pointpos,
                         boxpoints = boxpoints,
                         showlegend = FALSE,
                         marker = list(color = colors),
                         ...))
            } else{
              y = assay(object,assay_name)[row_name,]
              data <- purrr::pmap(list(split(y, groups), colors, levels(as.factor(groups))),
                                   function(j,k,l){
                                     list(y = j,
                                          x = l,
                                          name = ytitle,
                                          type = "box",
                                          pointpos = pointpos,
                                          boxpoints = boxpoints,
                                          marker = list(color = k),
                                          ...)
                                   })
            }
            
            # Make LocusSummary
            new("LocusSummary",
                data = data,
                layout = list(title = ytitle))
          })


setMethod("make_locus_summaries", c("SummarizedExperiment"),
          function(object,  
                   row_names,
                   assay_name = assayNames(object)[1],
                   ..., 
                   groups = NULL,
                   showlegend = !is.null(colors), 
                   colors = NULL,
                   boxpoints = c("all","Outliers","false"),
                   pointpos = 0,
                   ytitle = "Expression"){
            
            boxpoints = match.arg(boxpoints)
            
            if (is.null(colors)){
              colors <- "blue"
            }
            
            summaries <- purrr::map(seq_along(row_names), function(x){
              make_locus_summary(object, 
                                 row_names[x], 
                                 assay_name = assay_name,
                                 groups = groups, 
                                 showlegend = if (x == 1) showlegend else FALSE,
                                 legendgroup = "summary",
                                 colors = colors,
                                 boxpoints = boxpoints,
                                 pointpos = pointpos,
                                 ytitle = ytitle)
            })
            
            new("LocusSummaryList", as(summaries,"SimpleList"))
            
          })

#' make_summary_plotter
#' @param object SummarizedExperiment
#' @param assay_name name of assay to use
#' @param ... additional arguments
#' @param groups either vector of group assignments of name of column in object 
#' colData that corresponds to vector of group assignments
#' @param showlegend show the legend?
#' @param colors colors to use
#' @param boxpoints plot individual points?
#' @param pointpos relative position of points to boxes
#' @param ytitle name for yaxis
#' @export
setMethod("make_summary_plotter", c("SummarizedExperiment"),
          function(object,  
                   assay_name = assayNames(object)[1],
                   ..., 
                   groups = NULL,
                   showlegend = !is.null(colors), 
                   colors = NULL,
                   boxpoints = c("all","Outliers","false"),
                   pointpos = 0,
                   ytitle = "Expression"){
            
            boxpoints = match.arg(boxpoints)
            purrr::partial(make_locus_summaries,
                           object = object,
                           assay_name = assay_name,
                           ...,
                           groups = groups,
                           showlegend = showlegend,
                           boxpoints = boxpoints,
                           pointpos = pointpos,
                           ytitle = ytitle)
          })


setMethod(get_layout, "LocusSummary",
          function(object, yname, domain, anchor, ...){
            
            if (length(object@data) == 0) return(NULL)
            
            out <- list()
            
            # y axis settings
            out[[yname]] = modifyList(object@layout,
                                      list(zeroline = FALSE,
                                           domain = domain,
                                           anchor = gsub("xaxis","x",anchor),
                                           side = "right",
                                           ...))
            
            return(out)
            
            
          })

#' @export
setMethod(make_trace, signature = c(x = "LocusSummary"),
          definition = function(x, yax, xax, ...){
            lapply(x@data, function(y){ 
              y$yaxis = gsub("yaxis","y",yax)
              y$xaxis = gsub("xaxis","x",xax)
              y})
          })



