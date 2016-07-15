
#' single_coverage_plot
#' 
#' @param input either coverage matrix or input
#' @params regions the regions that will be plotted
#' @params region_names names of regions
#' @params binsize if inputs is filename, binsize for coverage
#' @params up if inputs is filename, basepairs upstream of tss to compute coverage for
#' @params down if inputs is filenames, basepairs downstream of tss to compute coverage for
#' @params row_agg optionally list of row aggregate profiles, otherwise computed
#' @params groups currently ignored
#' @params cvg_files currently accepts BigWig files
#' @params heatmap_normalization normalization to perform for heatmap, see Details
#' @params genome genome name, eg. "GRCm38" or "hg19"
#' @params org organism name, e.g. "mouse" or "human"
#' @return shiny app
#' @export
#' @import shiny
#' @author Alicia Schep
single_coverage_plot <- function(input, 
                                regions, 
                                region_names, 
                                binsize = 50,
                                up = 2000,
                                down = 2000,
                                row_agg = NULL,
                                groups = NULL,
                                cvg_files = NULL,
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

  if (!inherits(input,"matrix")){
    input <- make_coverage_matrix(inputs = input,
                                   ranges = regions,
                                   binsize = binsize,
                                   up = up,
                                   down = down)
  }
  
  input_normed = normalize_coverage_matrix(input, method = heatmap_normalization)
  
  if (is.null(row_agg)){
    #log10
    row_agg <- log10(rowSums(input) + 1)
  }
  

  single_coverage_shiny(input_normed, row_agg, cvg_files, regions, region_names, genome, org)
}

single_coverage_shiny <- function(M, signal, cvg_files, regions, region_names, genome, org){
  require(shiny)
  require(plotly)  
  
  # Check regions
  
  ui <- fluidPage(
    
    # Application title
    titlePanel("Heatmap Example"),
    
    fluidRow(
      column(4,
             fluidRow(
               selectInput("order",
                    "Order rows by",
                    choices = c("none","hclust","kmeans","signal"),
                    selected = "none"), 
               sliderInput("k",
                    "k = ",
                    min = 1,
                    max = 10,
                    value = 3,
                    step = 1)
             )),
      column(8, plotlyOutput("heat"))
    ),
    fluidRow(
      plotOutput("tracks")
    )
  )
  
  # Define server logic required to draw a histogram
  server <- function(input, output) {
    
    heat <-  reactive({single_coverage_heatmap(M,
                                     y = region_names, 
                                     row_order = input$order,
                                     k = input$k,
                                     signal = signal)})
    
    
    output$heat <- renderPlotly({
      heat()$plot()
    })
    
    output$tracks <- renderPlot({
      s <- event_data("plotly_click", source="HM")
      if(is.null(s) == T) return(NULL)
      ix <- heat()$row_order[s$y + 1]
      rang <- resize(regions[ix], width = 25000, fix = "center")
      de <- igisExonlist(rang, genome = genome, org = org)
      grl <- getCoverageInRange(cvg_files, rang, names = names(cvg_files))
      plotGeneCoverage(grl, de, rang, genome = genome,symbol=region_names[ix])
    })
    
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server, options = list(height = 1400))
  
}


