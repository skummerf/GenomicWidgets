
HEATMAP_SOURCE <- "HM"

#' heatmap_click
#' 
#' Function to bind heatmap to click event
#' @param heatmap heatmap object
#' @param x list-like object corresponding to rows of heatmap
#' @return function that returns appropriate element of x based on row of 
#' heatmap clicked within \code{\link[shiny]{shinyApp}}
#' @return returns function
#' @author Alicia Schep 
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
#' Function for making shiny app linking a heatmap to genome tracks
#' 
#' @param heatmap IHeatmap object
#' @param track_params output from \code{\link{set_track_parameters}}
#' @param link link function, linking rows of heatmap to input to track function
#' generator, result from \code{\link{heatmap_click}}
#' @param title Title of shiny app
#' @param options to pass to shiny
#' @return Shiny application
#' @author Alicia Schep and Justin Finkle
#' 
#' @export
#' 
heatmap_to_tracks_shiny <- function(heatmap, 
                                    track_params,
                                    link,
                                    title = "Heatmap linked to Genome Tracks",
                                    options = list(height = 1400)){
  
  if (!(requireNamespace("shiny"))) stop("Must have shiny package installed!")
  
  # Check regions
  
  ui <- shiny::fluidPage(
    
    # Application title
    shiny::titlePanel(title),
    
    shiny::fluidRow(
      iheatmaprOutput("heat")
    ),
    shiny::fluidRow(
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
      if (is.null(linker)) NULL else plot_tracks(linker, track_params) 
    })
    
  }
  
  # Run the application 
  shiny::shinyApp(ui = ui, server = server, options = options)
  
}

