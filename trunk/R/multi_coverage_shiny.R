
multi_coverage_plot <- function(inputs, 
                                regions, 
                                region_names, 
                                binsize = 50,
                                up = 2000,
                                down = 2000,
                                col_aggs = NULL,
                                row_aggs = NULL,
                                groups = NULL,
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
                                ...){
  
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
    row_aggs <- lapply(inputs, function(x) log10(rowSums(x)))
  }
                       
  multi_coverage_shiny(inputs_normed, 
                       colnames(inputs[[1]]),
                       region_names, 
                       col_aggs = col_aggs, 
                       row_aggs = row_aggs,
                       ranges = regions,
                       groups = groups)
  
}




#Probably want a wrapper around this one
multi_coverage_shiny <- function(mats, 
                                 positions, 
                                 region_names, 
                                 col_aggs,
                                 row_aggs,
                                 groups = NULL, 
                                 ranges = NULL){
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
        
      )#,
      
    #fluidRow( plotOutput("tracks"))
    
    )
  
  # Define server logic required to draw a histogram
  server <- function(input, output) {
    
    output$agg <- renderPlotly({
       
      plot_ly(x = rep(positions, length(mats)),
              y = unlist(col_aggs, use.names = FALSE),
              color = rep(names(mats), each = length(positions)),
              text = rep(names(mats), each = length(positions)),
              mode = "lines",
              source = "agg") %>% layout(hovermode = "closest",
                                         xaxis = list(title = "Position"),
                                         yaxis = list(title = "ChIP Signal"),
                                         showlegend = TRUE)
      
    })
    
    output$heat <- renderPlotly({
      s <- event_data("plotly_click", source="agg")
      if (length(s) != 0){
        sel <- s$curve + 1
        sig <- row_aggs[[sel]]
        highsig <- which(rank(sig) > length(sig) - 1000)
        single_coverage_heatmap(mats[[sel]][highsig,], x = positions, y = region_names, row_order = "signal", 
                                signal = sig[highsig]) %>% layout(title = names(mats)[sel])
        
      } else{
        plotly_empty()
      }

    })
    # 
    # output$track <- renderPrint({
    #   s1 <- event_data("plotly_click", source="agg")
    #   if (length(s1) != 0){
    #     s <- event_data("plotly_click", source="HM")
    #     if (length(s) != 0){
    #       str(s2$y)
    #     }
    #   } else{
    #     ""
    #   }
    # })
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server)
  
}