
HEATMAP_SOURCE = "HM"

#' heatmap_click
#' 
#' @export
heatmap_click <- function(heatmap, ranges, width = 50000){
  ranges <- resize(ranges, fix = "center", width = width)
  out <- function(){
    s <- event_data("plotly_click", source = HEATMAP_SOURCE)
    if(is.null(s)) return(NULL)
    yaxis <- heatmap$plots %>% filter(pid == s$curve) %>% select(yid) %>% unlist()
    if (yaxis != "y") return(NULL)
    ro <- heatmap$yaxes %>% filter(yid == yaxis) %>% select(row_order) %>% unlist()
    ix <- ro[s$y + 1]
    return(ranges[ix])
  }
  return(out)
}

#' heatmap_to_tracks_shiny
#' 
#' @export
heatmap_to_tracks_shiny <- function(heatmap, 
                                    track_function,
                                    link,
                                    title = "Heatmap linked to Genome Tracks",
                                    options = list(height = 1400)){
  require(shiny)
  require(plotly)  
  
  # Check regions
  
  ui <- fluidPage(
    
    # Application title
    titlePanel(title),
    
    fluidRow(
      plotlyOutput("heat")
    ),
    fluidRow(
      plotOutput("tracks")
    )
  )
  
  # Define server logic required to draw a histogram
  server <- function(input, output) {
    
    output$heat <- renderPlotly({
      heatmap %>% as_plotly(source = HEATMAP_SOURCE)
    })
    
    output$tracks <- renderPlot({
      linker <- link() 
      track_function(linker)  
    })
    
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server, options = options)
  
}

#' heatmap_to_browserly_shiny
#' 
#' @export
heatmap_to_browserly_shiny <- function(heatmap, 
                                    track_function,
                                    link,
                                    title = "Heatmap linked to Genome Tracks",
                                    options = list(height = 1400)){
  require(shiny)
  require(plotly)  
  
  # Check regions
  
  ui <- fluidPage(
    
    # Application title
    titlePanel(title),
    
    fluidRow(
      plotlyOutput("heat")
    ),
    fluidRow(
      plotlyOutput("tracks")
    )
  )
  
  # Define server logic required to draw a histogram
  server <- function(input, output) {
    
    output$heat <- renderPlotly({
      heatmap %>% as_plotly(source = HEATMAP_SOURCE)
    })
    
    output$tracks <- renderPlotly({
      linker <- link() 
      if (is.null(linker)) return(NULL)
      track_function(linker)  
    })
    
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server, options = options)
  
}
