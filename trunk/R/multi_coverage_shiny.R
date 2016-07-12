
multi_coverage_plot <- function(inputs, 
                                regions, 
                                region_names, 
                                binsize = 50,
                                up = 2000,
                                down = 2000,
                                col_aggs = NULL,
                                row_aggs = NULL,
                                groups = NULL,
                                cvg_files = NULL,
                                aggregate_normalization = c("localRms", 
                                                            "localMean", 
                                                            "localNonZeroMean", 
                                                            "PercentileMax", 
                                                            "scalar", 
                                                            "none"),
                                heatmap_normalization = c("localRms", 
                                                          "localMean", 
                                                          "localNonZeroMean", 
                                                          "PercentileMax", 
                                                          "scalar", 
                                                          "none"),
                                genome = "GRCm38",
                                org = "mouse",
                                ...){
  
  heatmap_normalization = match.arg(heatmap_normalization)
  aggregate_normalization = match.arg(aggregate_normalization)
  
  if (!inherits(inputs[[1]],"matrix")){
    inputs <- make_coverage_matrix(inputs = inputs,
                       ranges = regions,
                       binsize = binsize,
                       up = up,
                       down = down)
  }
  
  if (is.null(col_aggs)){
    inputs_normed = normalize_coverage_matrix(inputs, method = aggregate_normalization)
    col_aggs = lapply(inputs_normed, colSums)
  }
  
  if (heatmap_normalization != aggregate_normalization){
    inputs_normed = normalize_coverage_matrix(inputs, method = heatmap_normalization)
  }
  
  if (is.null(row_aggs)){
    #log10
    row_aggs <- lapply(inputs, function(x) log10(rowSums(x) + 1))
  }
  
  heatmaps <- BiocParallel::bplapply(1:length(inputs), function(x){
    sig <- row_aggs[[x]]
    if (length(sig) > 500){
      highsig <- which(sig > quantile(sig,(length(sig) - 500)/length(sig)))
    } else{
      highsig <- 1:length(sig)
    }
    out <- single_coverage_heatmap(inputs_normed[[x]][highsig,],
                                   y = region_names[highsig], 
                                   row_order = "signal", 
                            signal = sig[highsig]) 
    out <- c(out, list(row_mapping = highsig[out$row_order]))
    return(out)
  })
                       
  multi_coverage_shiny(heatmaps,
                       positions = colnames(inputs[[1]]),
                       region_names,
                       col_aggs = col_aggs,
                       regions = regions,
                       cvg_files = cvg_files,
                       genome = genome,
                       org = org)
  
}




#Probably want a wrapper around this one
multi_coverage_shiny <- function(heatmaps, 
                                 positions,
                                 region_names, 
                                 col_aggs,
                                 cvg_files = NULL,
                                 regions = NULL,
                                 genome = NULL,
                                 org = NULL){
  require(shiny)
  require(plotly)  
  
  #aggs <- lapply(mats, colMeans)
  
  # Check regions
  
  ui <- fluidPage(
    
    # Application title
    titlePanel("TSS multiple samples"),
    
    #
     fluidRow(
        column(6, plotlyOutput("agg")),
        column(6, plotlyOutput("heat"))
        
      ),
      
    fluidRow( column(12,plotOutput("tracks")))
    
    )
  
  # Define server logic required to draw a histogram
  server <- function(input, output) {
    
    output$agg <- renderPlotly({
       
      plot_ly(x = rep(positions, length(col_aggs)),
              y = unlist(col_aggs, use.names = FALSE),
              color = rep(names(col_aggs), each = length(positions)),
              text = rep(names(col_aggs), each = length(positions)),
              mode = "lines",
              source = "agg") %>% layout(hovermode = "closest",
                                         xaxis = list(title = "Position"),
                                         yaxis = list(title = "ChIP Signal"),
                                         showlegend = TRUE)
      
    })
    
    
    
    output$heat <- renderPlotly({
      s1 <- event_data("plotly_click", source="agg")
      if(is.null(s1) == T) return(NULL)
      
      sel <- s1$curve + 1
      heatmaps[[sel]]$plot() %>% layout(title = names(col_aggs)[sel])
      

    })

    output$tracks <- renderPlot({
      s1 <- event_data("plotly_click", source="agg")
      if(is.null(s1) == T) return(NULL)
      sel <- s1$curve + 1
      s2 <- event_data("plotly_click", source="HM")
        
        # If NULL dont do anything
      if(is.null(s2) == T) return(NULL)
        
      ix <- heatmaps[[sel]]$row_mapping[s2$y + 1]
      rang <- resize(regions[ix], width = 50000, fix = "center")
      de <- igisExonlist(rang, genome = genome, org = org)
      grl <- getCoverageInRange(cvg_files, rang, names = names(cvg_files))
      plotGeneCoverage(grl, de, rang, genome = genome,symbol=region_names[ix])

    })
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server, options = list(height = 1200))
  
}