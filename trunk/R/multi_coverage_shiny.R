

#Probably want a wrapper around this one

multi_coverage_shiny <- function(mats, positions){
  require(shiny)
  require(plotly)  
  
  aggs <- lapply(mats, colMeans)
  
  # Check regions
  
  ui <- fluidPage(
    
    # Application title
    titlePanel("TSS multiple samples"),
    
      # Show a plot of the generated distribution
      fluidRow(
        column(6, plotlyOutput("agg")),
        column(6, plotlyOutput("heat"))
        
      )
    )
  
  # Define server logic required to draw a histogram
  server <- function(input, output) {
    
    output$agg <- renderPlotly({
       
      plot_ly(x = rep(positions, length(mats)),
              y = unlist(aggs, use.names = FALSE),
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
        single_coverage_heatmap(mats[[sel]], x = positions, row_order = "signal", 
                                signal = rowMeans(mats[[sel]])) %>% layout(title = names(mats)[sel])
        
      } else{
        plotly_empty()
      }

    })
    
    
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server)
  
}