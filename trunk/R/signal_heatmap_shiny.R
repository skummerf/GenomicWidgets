#' rna_chip_tracks_plot
#' 
#' @param mat matrix
#' @param x x axis labels
#' @param y y axis labels
#' @param row_order row order method
#' @param col_order col order method
#' @param row_k k to use for kmeans clustering or cutting heirarchical clustering
#' @param col_k k to use for kmeans clustering or cutting heirarchical clustering
#' @param row_groups pre-determined groups for rows
#' @param col_groups pre-determined groups for col
#' @param row_clust_dist distance function to uses for row clustering
#' @param col_clust_dist distance function to uses for col clustering
#' @param name name of colorbar
#' @params regions the regions that will be plotted
#' @params region_names names of regions
#' @params cvg_files currently accepts BigWig files
#' @params genome genome name, eg. "GRCm38" or "hg19"
#' @params org organism name, e.g. "mouse" or "human"
#' @return shiny app
#' @export
#' @import shiny
#' @author Alicia Schep
rna_chip_tracks_plot <- function(mat, 
                                 x = chipVis:::default_x(mat),
                                 y = chipVis:::default_y(mat),                   
                                 row_order = c("none","hclust","kmeans","groups"),
                                 col_order = c("none","hclust","kmeans","groups"),
                                 row_groups = NULL,
                                 col_groups = NULL,
                                 row_anno = NULL,
                                 col_anno = NULL,
                                 row_k = NULL,
                                 col_k = NULL,
                                 row_clust_dist = stats::dist,
                                 col_clust_dist = stats::dist,
                                 name = "Row Z-score",
                                 scale = TRUE,
                                 cvg_files, 
                                 regions, 
                                 region_names, 
                                 genome, 
                                 org){
  require(shiny)
  require(plotly)  
  
  # Check regions
  
  ui <- fluidPage(
    
    # Application title
    titlePanel("Heatmap Example"),
    
    fluidRow(
      column(4,
             fluidRow(
               selectInput("row_order",
                           "Order rows by",
                           choices = c("none","hclust","kmeans","signal"),
                           selected = "none"), 
               sliderInput("row_k",
                           "k for row clustering",
                           min = 1,
                           max = 10,
                           value = 3,
                           step = 1),
               selectInput("col_order",
                           "Order columns by",
                           choices = c("none","hclust","kmeans","signal"),
                           selected = "none"), 
               sliderInput("col_k",
                           "k for column clustering",
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
    
    heat <-  reactive({cluster_heatmap(mat, 
                                       x = x,
                                       y = y,
                                       row_order = input$row_order,
                                       row_k = input$row_k,
                                       col_order = input$col_order,
                                       col_k = input$col_k,
                                       row_anno = row_anno,
                                       col_anno = col_anno,
                                       row_clust_dist = row_clust_dist,
                                       col_clust_dist = col_clust_dist,
                                       name = name,
                                       source = "HM",
                                       scale = scale)})
    
    output$heat <- renderPlotly({
      heat()$plot() %>% layout(margin = list(b= 100))
    })
    
    output$tracks <- renderPlot({
      s <- event_data("plotly_click", source="HM")
      if(is.null(s) == T) return(NULL)
      ix <- heat()$row_order[s$y + 1]
      rang <- GenomicRanges::resize(regions[ix], width = 25000, fix = "center")
      de <- igisExonlist(rang, genome = genome, org = org)
      grl <- getCoverageInRange(cvg_files, rang, names = names(cvg_files))
      plotGeneCoverage(grl, de, rang, genome = genome,symbol=region_names[ix])
    })
    
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server, options = list(height = 1400))
  
}
