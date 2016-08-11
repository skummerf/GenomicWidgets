
HEATMAP_SOURCE = "HM"

heatmap_click <- function(heatmap, ranges, width = 50000){
  ranges <- resize(ranges, fix = "center", width = width)
  out <- function(){
    s <- event_data("plotly_click", source = HEATMAP_SOURCE)
    if(is.null(s)) return(NULL)
    if (heatmap$curve_map[[s$curve + 1]]$yaxis != "yaxis") return(NULL)
    ix <- heatmap$row_order[s$y + 1]
    return(ranges[ix])
  }
  return(out)
}

tracks_placeholder <- function(cvg_files, txdb, tx_data, genome, org, scaling_factor, ...){
  out <- function(range){
    if (is.null(range)) return(NULL)
    de <- get_tx_annotation(db_object = txdb, 
                            range = range, 
                            tx_data = tx_data,
                            no_introns=TRUE)
    grl <- get_coverage_in_range(bwList = cvg_files, 
                                 target_range = range, 
                                 names = names(cvg_files), 
                                 cvg_scaling  = scaling_factor)
    plot_track_view(cvg_list =  grl, 
                    target_range = range,
                    exon_data = de,
                    genome = genome, ...)
  }
  return(out)
}

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
      heatmap %>% plot_iHeatmap(source = HEATMAP_SOURCE)
    })
    
    output$tracks <- renderPlot({
      linker <- link() 
      track_function(linker)  
    })
    
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server, options = options)
  
}
