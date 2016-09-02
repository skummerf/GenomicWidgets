AGGREGATE_SOURCE = "Agg"

#' aggregate_shiny
#' 
#' @param sample_groups vector of groups to which matrices belong
#' @param colors either RColorbrewer pallete name or vector of colors
#' @param xlab axis name for x axis
#' @param pct percentile argument if PercentileMax chosen as scale_method
#' @param scale_factors scale_factors if "scalar" chosen as scale_method
#' @param mats 
#' @param aggregate_scale_method 
#' @param showlegend show the legend?
#' @param x 
#' @param y 
#' @param include 
#' @param include_method 
#' @param row_order 
#' @param k 
#' @param clust_dist 
#' @param signal 
#' @param plot_signal 
#' @param name 
#' @param signal_name 
#' @param summary 
#' @param heatmap_scale_method 
#' @param ticktext 
#' @param start 
#' @param end 
#' @param font 
#' @param title 
#' @param options 
#' @param ... 
#'
#' @details  More details will go here.
#' @export
#' @author Alicia Schep
aggregate_shiny <- function(mats,
                            sample_groups = NULL,
                            colors = "Dark2",
                            aggregate_scale_method = c("none",
                                             "localRms", 
                                             "localMean", 
                                             "localNonZeroMean", 
                                             "PercentileMax", 
                                             "scalar"),
                            showlegend = TRUE,
                            x = iHeatmapR:::default_x(mats[[1]]),
                            y = iHeatmapR:::default_y(mats[[1]]),
                            include = 1000,
                            include_method = c("signal","first","random"),
                            row_order = c("none","hclust","kmeans","groups","signal"),
                            k = NULL,
                            row_groups = NULL,
                            clust_dist = stats::dist,
                            signal = lapply(mats, rowMeans, na.rm = TRUE),
                            plot_signal = TRUE,
                            name = "Coverage",
                            signal_name = "Avg. Coverage",
                            summary = TRUE,
                            heatmap_scale_method = c("localRms", 
                                             "localMean", 
                                             "localNonZeroMean", 
                                             "PercentileMax", 
                                             "scalar", 
                                             "none"),
                            pct = 0.95,
                            scale_factors = rep(1, length(mats)), 
                            ticktext = TRUE,
                            start = x[1],
                            end = default_end(x),     
                            xlab = "Position",
                            font = list(size = 8),
                            title = "Aggregate Plot Linked to Heatmap",
                            options = list(height = 1400),
                            ...){
  
  
  # Add Argument Check
  
  heatmap_scale_method = match.arg(heatmap_scale_method)
  aggregate_scale_method = match.arg(aggregate_scale_method)
  include_method = match.arg(include_method)
  row_order = match.arg(row_order)
  
  agg <- aggregate_profile_plot(mats, 
                                groups = sample_groups, 
                                positions = as.numeric(x), 
                                ylab = name,
                                xlab = xlab, 
                                scale_method = aggregate_scale_method,
                                pct = pct,
                                scale_factors = scale_factors, 
                                showlegend = showlegend, 
                                colors = colors,
                                source = AGGREGATE_SOURCE)
  

  heatmaps <- lapply(seq_along(mats), function(ix){
    single_coverage_heatmap(mats[[ix]], 
                                x = x, 
                                y = y, 
                                include = include, 
                                include_method = include_method, 
                                row_order = row_order, 
                                k = k, 
                                groups = row_groups, 
                                clust_dist = clust_dist, 
                                name = name, 
                                signal_name = signal_name, 
                                scale_method = heatmap_scale_method, 
                                summary = summary, 
                                source = "HM", 
                                plot_signal = plot_signal, 
                                scale_factor = scale_factors[ix], 
                                pct = pct, 
                                signal = signal[[ix]],
                                ticktext = ticktext, 
                                start = start, 
                                end = end, 
                                xlab = xlab, 
                                font = font,
                                ...) 
  })
  
  aggregate_to_heatmap_shiny(agg, heatmaps, title, options)
  
}

#' aggregate_shiny_plus
#' 
#' @param sample_groups vector of groups to which matrices belong
#' @param colors either RColorbrewer pallete name or vector of colors
#' @param xlab axis name for x axis
#' @param pct percentile argument if PercentileMax chosen as scale_method
#' @param scale_factors scale_factors if "scalar" chosen as scale_method
#' @param mats 
#' @param aggregate_scale_method 
#' @param showlegend show the legend?
#' @param x 
#' @param y 
#' @param include 
#' @param include_method 
#' @param row_order 
#' @param k 
#' @param clust_dist 
#' @param signal 
#' @param plot_signal 
#' @param name 
#' @param signal_name 
#' @param summary 
#' @param heatmap_scale_method 
#' @param ticktext 
#' @param start 
#' @param end 
#' @param font 
#' @param title 
#' @param options 
#' @param ... 
#'
#' @details  More details will go here.
#' @export
#' @author Alicia Schep
aggregate_shiny_plus <- function(mats,
                                 regions,
                                 track_function,
                                 region_widths = 50000,
                            sample_groups = NULL,
                            colors = "Dark2",
                            aggregate_scale_method = c("none",
                                                       "localRms", 
                                                       "localMean", 
                                                       "localNonZeroMean", 
                                                       "PercentileMax", 
                                                       "scalar"),
                            showlegend = TRUE,
                            x = iHeatmapR:::default_x(mats[[1]]),
                            y = iHeatmapR:::default_y(mats[[1]]),
                            include = 1000,
                            include_method = c("signal","first","random"),
                            row_order = c("none","hclust","kmeans","groups","signal"),
                            k = NULL,
                            row_groups = NULL,
                            clust_dist = stats::dist,
                            signal = lapply(mats, rowMeans, na.rm = TRUE),
                            plot_signal = TRUE,
                            name = "Coverage",
                            signal_name = "Avg. Coverage",
                            summary = TRUE,
                            heatmap_scale_method = c("localRms", 
                                                     "localMean", 
                                                     "localNonZeroMean", 
                                                     "PercentileMax", 
                                                     "scalar", 
                                                     "none"),
                            pct = 0.95,
                            scale_factors = rep(1, length(mats)), 
                            ticktext = TRUE,
                            start = x[1],
                            end = default_end(x),     
                            xlab = "Position",
                            font = list(size = 8),
                            title = "Aggregate Plot Linked to Heatmap",
                            options = list(height = 1400),
                            ...){
  
  
  # Add Argument Check
  
  heatmap_scale_method = match.arg(heatmap_scale_method)
  aggregate_scale_method = match.arg(aggregate_scale_method)
  include_method = match.arg(include_method)
  row_order = match.arg(row_order)
  
  agg <- aggregate_profile_plot(mats, 
                                groups = sample_groups, 
                                positions = as.numeric(x), 
                                ylab = name,
                                xlab = xlab, 
                                scale_method = aggregate_scale_method,
                                pct = pct,
                                scale_factors = scale_factors, 
                                showlegend = showlegend, 
                                colors = colors,
                                source = AGGREGATE_SOURCE)
  
  
  heatmaps <- lapply(seq_along(mats), function(ix){
    single_coverage_heatmap(mats[[ix]], 
                            x = x, 
                            y = y, 
                            include = include, 
                            include_method = include_method, 
                            row_order = row_order, 
                            k = k, 
                            groups = row_groups, 
                            clust_dist = clust_dist, 
                            name = name, 
                            signal_name = signal_name, 
                            scale_method = heatmap_scale_method, 
                            summary = summary, 
                            source = "HM", 
                            plot_signal = plot_signal, 
                            scale_factor = scale_factors[ix], 
                            pct = pct, 
                            signal = signal[[ix]],
                            ticktext = ticktext, 
                            start = start, 
                            end = end, 
                            xlab = xlab, 
                            font = font,
                            ...) 
  })
  
  links <- lapply(heatmaps, heatmap_click, regions, width = region_widths)
  
  if (inherits(track_function,"browserly")){
    out <- aggregate_to_heatmap_to_browserly_shiny(agg, 
                                          heatmaps, 
                                          track_function,
                                          links,
                                          title, 
                                          options)
  } else{
    out <- aggregate_to_heatmap_to_tracks_shiny(agg, 
                                                   heatmaps, 
                                                   track_function,
                                                   links,
                                                   title, 
                                                   options)
  }
  return(out)
  
}


#' aggregate_to_heatmap_shiny
#' 
#' @export
aggregate_to_heatmap_shiny <- function(aggregate_plot,
                                       heatmaps, 
                                        title = "Aggregate Plot Linked to Heatmap",
                                        options = list(height = 1400)){
  require(shiny)
  require(plotly)  
  
  # Check regions
  
  ui <- fluidPage(
    
    # Application title
    titlePanel(title),
    
    fluidRow(
      column(6, plotlyOutput("agg")),
      column(6, plotlyOutput("heat"))
    )
  )
  
  # Define server logic required to draw a histogram
  server <- function(input, output) {
    
    output$agg <- renderPlotly({
      aggregate_plot
    })
    
    output$heat <- renderPlotly({
      s <- event_data("plotly_click", source = AGGREGATE_SOURCE)
      if(is.null(s)) return(NULL)
      sel <- s$curve + 1
      heatmaps[[sel]] %>% as_plotly()
    })
    
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server, options = options)
  
}

#' heatmap_to_browserly_shiny
#' 
#' @export
aggregate_to_heatmap_to_browserly_shiny <- function(aggregate_plot,
                                                    heatmaps, 
                                       track_function,
                                       link_functions,
                                       title = "Heatmap linked to Genome Tracks",
                                       options = list(height = 1400)){
  require(shiny)
  require(plotly)  
  
  # Check regions
  
  ui <- fluidPage(
    
    # Application title
    titlePanel(title),
    
    fluidRow(
      column(6, plotlyOutput("agg")),
      column(6, plotlyOutput("heat"))
    ),
    fluidRow(
      plotlyOutput("tracks")
    )
  )
  
  # Define server logic required to draw a histogram
  server <- function(input, output) {
    
    output$agg <- renderPlotly({
      aggregate_plot
    })
    
    
    agg_linker <- function()({
      s <- event_data("plotly_click", source = AGGREGATE_SOURCE)
      if(is.null(s)){
        return(NULL)
      } 
      else{
        return(s$curve + 1)
      }
    })
    
    #sel <- reactive({})
    
    output$heat <- renderPlotly({
      sel <- agg_linker()
      if (is.null(sel)) return(NULL)
      heatmaps[[sel]] %>% as_plotly(source = HEATMAP_SOURCE)
    })
    
    output$tracks <- renderPlotly({
      sel <- agg_linker()
      if (is.null(sel)) return(NULL)
      linker <- link_functions[[sel]]() 
      if (is.null(linker)) return(NULL)
      track_function(linker)  
    })
    
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server, options = options)
  
}

#' heatmap_to_browserly_shiny
#' 
#' @export
aggregate_to_heatmap_to_tracks_shiny <- function(aggregate_plot,
                                                    heatmaps, 
                                                    track_function,
                                                    link_functions,
                                                    title = "Heatmap linked to Genome Tracks",
                                                    options = list(height = 1400)){
  require(shiny)
  require(plotly)  
  
  # Check regions
  
  ui <- fluidPage(
    
    # Application title
    titlePanel(title),
    
    fluidRow(
      column(6, plotlyOutput("agg")),
      column(6, plotlyOutput("heat"))
    ),
    fluidRow(
      plotOutput("tracks")
    )
  )
  
  # Define server logic required to draw a histogram
  server <- function(input, output) {
    
    output$agg <- renderPlotly({
      aggregate_plot
    })
    
    
    agg_linker <- function()({
      s <- event_data("plotly_click", source = AGGREGATE_SOURCE)
      if(is.null(s)){
        return(NULL)
      } 
      else{
        return(s$curve + 1)
      }
    })
    
    #sel <- reactive({})
    
    output$heat <- renderPlotly({
      sel <- agg_linker()
      if (is.null(sel)) return(NULL)
      heatmaps[[sel]] %>% as_plotly(source = HEATMAP_SOURCE)
    })
    
    output$tracks <- renderPlot({
      sel <- agg_linker()
      if (is.null(sel)) return(NULL)
      linker <- link_functions[[sel]]() 
      if (is.null(linker)) return(NULL)
      track_function(linker)  
    })
    
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server, options = options)
  
}



