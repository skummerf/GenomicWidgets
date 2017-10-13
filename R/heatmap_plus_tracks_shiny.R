
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
#' @examples 
#' 
#' library(GenomicRanges)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' 
#'   
#' ## We'll also read in some track data to plot
#' genomation_dir <- system.file("extdata", package = "genomationData")
#' samp.file <- file.path(genomation_dir,'SamplesInfo.txt')
#' samp.info <- read.table(samp.file, header=TRUE, sep="\t", 
#'                         stringsAsFactors = FALSE)
#' samp.info$fileName <- file.path(genomation_dir, samp.info$fileName)
#' ctcf.peaks = genomation::readBroadPeak(system.file("extdata",
#'                          "wgEncodeBroadHistoneH1hescCtcfStdPk.broadPeak.gz",
#'                           package = "genomationData"))
#' ctcf.peaks = ctcf.peaks[seqnames(ctcf.peaks) == "chr21"]
#' 
#' ## resize peaks to size 1000
#' ctcf.peaks = resize(ctcf.peaks, width = 10000, fix = "center")
#' 
#' ## Make track plotter using summary parametrs
#' 
#' track_params <- set_track_parameters(samp.info$fileName[1:3], 
#'   annotation = TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'   track_names = samp.info$sampleName[1:3], 
#'   share_y = TRUE)
#' 
#' # Make coverage heamap
#' 
#' ctcf_mats <- make_coverage_matrix(samp.info$fileName[1:5], 
#'                                   ctcf.peaks, 
#'                                   input_names = samp.info$sampleName[1:5],
#'                                   up = 250, 
#'                                   down = 250, 
#'                                   binsize = 25)
#'                                   
#' hm <- coverage_heatmap(ctcf_mats, "Ctcf") 
#' 
#' link_fn <- heatmap_click(hm, ctcf.peaks)
#'       
#' if (interactive()){
#'   heatmap_to_tracks_shiny(hm, track_params, link_fn)
#' }   
#' 
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
#' @examples 
#' 
#' library(GenomicRanges)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' 
#'   
#' ## We'll also read in some track data to plot
#' genomation_dir <- system.file("extdata", package = "genomationData")
#' samp.file <- file.path(genomation_dir,'SamplesInfo.txt')
#' samp.info <- read.table(samp.file, header=TRUE, sep="\t", 
#'                         stringsAsFactors = FALSE)
#' samp.info$fileName <- file.path(genomation_dir, samp.info$fileName)
#' ctcf.peaks = genomation::readBroadPeak(system.file("extdata",
#'                          "wgEncodeBroadHistoneH1hescCtcfStdPk.broadPeak.gz",
#'                           package = "genomationData"))
#' ctcf.peaks = ctcf.peaks[seqnames(ctcf.peaks) == "chr21"]
#' 
#' ## resize peaks to size 1000
#' ctcf.peaks = resize(ctcf.peaks, width = 10000, fix = "center")
#' 
#' ## Make track plotter using summary parametrs
#' 
#' track_params <- set_track_parameters(samp.info$fileName[1:3], 
#'   annotation = TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'   track_names = samp.info$sampleName[1:3], 
#'   share_y = TRUE)
#' 
#' # Make coverage heamap
#' 
#' ctcf_mats <- make_coverage_matrix(samp.info$fileName[1:5], 
#'                                   ctcf.peaks, 
#'                                   input_names = samp.info$sampleName[1:5],
#'                                   up = 250, 
#'                                   down = 250, 
#'                                   binsize = 25)
#'                                   
#' hm <- coverage_heatmap(ctcf_mats, "Ctcf") 
#' 
#' link_fn <- heatmap_click(hm, ctcf.peaks)
#'       
#' if (interactive()){
#'   heatmap_to_tracks_shiny(hm, track_params, link_fn)
#' }   
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

