
HEATMAP_SOURCE = "HM"

#' heatmap_click
#' 
#' @export
heatmap_click <- function(heatmap, x){
  out <- function(){
    s <- iheatmapr_event(heatmap, "click")
    if(is.null(s)) return(NULL)
    ix <- s$row
    return(x[ix])
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
  requireNamespace(shiny)
  
  # Check regions
  
  ui <- fluidPage(
    
    # Application title
    titlePanel(title),
    
    fluidRow(
      iheatmaprOutput("heat")
    ),
    fluidRow(
      GenomicWidgetsOutput("tracks")
    )
  )
  
  # Define server logic required to draw a histogram
  server <- function(input, output) {
    
    output$heat <- renderIheatmap({
      heatmap
    })
    
    output$tracks <- renderGenomicWidgets({
      linker <- link() 
      if (is.null(linker)) NULL else track_function(linker) 
    })
    
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server, options = options)
  
}

