#' Title
#'
#' @param cvg_files 
#' @param genome 
#' @param tx_data 
#' @param hm_thresh 
#' @param type 
#' @param cvg_scaling 
#'
#' @return
#' @export
#'
#' @examples
make_browserly_function <- function(cvg_files, 
                                    tx_data,
                                    sample_names = names(cvg_files),
                                    hm_thresh = 4,
                                    type = NULL,
                                    cvg_scaling = rep(1, length(cvg_files))){
  plot_browserly <- function(target_range){
    cvg_gr <- make_coverage_tracks(inputs = cvg_files,
                                      target_range = target_range, 
                                      sample_names = names(cvg_files), 
                                      scaling_factors = cvg_scaling)
    b_plot <- plot_single_locus(target_range = target_range, 
                                tx_data = tx_data,
                                cvg = cvg_gr,
                                hm_thresh = hm_thresh,
                                type = type)
    return(b_plot)
  }
  class(plot_browserly) <- c(class("plot_browserly"), "browserly")
  return(plot_browserly)
}


#' Title
#'
#' @param target_range 
#' @param tx_data 
#' @param cvg_L 
#' @param hm_thresh 
#' @param stacking 
#' @param sync_y 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot_multiple_genes <- function(genes,
                                centered_cvg,
                                centered_tx,
                                hm_thresh=4,
                                stacking = 'dense',
                                type = 'scatter',
                                fill = 'none',
                                mode = 'lines',
                                sync_y = FALSE,
                                xlabel = 'TSS',
                                ...){
  # Make transcript tracks
  tx_tracks <- lapply(genes, function(gene){
    tx_info <- centered_tx[[gene]]
    tx_info <- set_tx_level(tx_info, stacking = stacking)
    tx_track <- browserly_tx_track(tx_info = tx_info, 
                                   track_name = paste0(gene,'_Annotation'))
    return(tx_track[[1]])
  })
  
  cvg_tracks <- make_subplots(grl = cvg_list, 
                              type = type, 
                              fill = fill, 
                              mode = mode,
                              legend = 'first',
                              ...)
  
  # Order tracks, alternating coverage then annotation
  plot_order <- order(rep(seq_along(genes), 2))
  plots <- c(cvg_tracks, tx_tracks)[plot_order]
  sp <- subplot(plots, shareX = TRUE, nrows=length(plots))
  sp_info <- get_subplot_ax_info(sp)
  
  display_range <- c(min(min(start(centered_cvg))), max(max(end(centered_cvg))))
  
  for(gene in genes){
    # Get the axis on which annotations are made
    tx_info <- tx_list[[gene]]
    tx_info <- set_tx_level(tx_info, stacking = 'dense')
    grt_ax <- filter(sp_info, grepl(gene, subplot_name), 
                     is_annotation == TRUE)$yaxis[[1]]
    # Add the annotation shapes
    sp <- add_tx_shapes(plotly_obj = sp, 
                        tx_info = tx_info, 
                        target_range = range(tx_info), 
                        grt_ax = grt_ax)
    
  }
  # Adjust the layout parameters
  heights <- rep(c((0.25/length(genes)), (0.75/length(genes))), length(genes))
  sp <- adjust_y_domains(plotly_obj = sp,
                         ax_info = sp_info,
                         heights = heights)
  sp <- modify_y(plotly_obj = sp, 
                 ax_info = sp_info, 
                 grt_ax = grt_ax, 
                 type = type, 
                 sync_y = sync_y, 
                 native_title = TRUE)
  
  # Minor Layout tweaks
  sp$x$layout$xaxis$range <- display_range
  sp$x$layout$xaxis$title <- xlabel
  sp$x$layout$xaxis$zeroline <- FALSE
  sp$x$layout$margin <- list(t=0, b = 30, r = 0)
  sp$x$layout$hovermode <- 'compare'
  
  return(sp)
}

#' Title
#'
#' @param tx_data 
#' @param cvg
#' @param hm_thresh 
#' @param type 
#' @param target_range 
#' @param stacking 
#' @param sync_y 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot_single_locus <- function(target_range, 
                              tx_data, 
                              cvg,
                              hm_thresh = 4,
                              stacking = 'dense',
                              sync_y = TRUE,
                              type = NULL,
                              xlabel = NULL,
                              ...){
  # Make GeneRegion plot
  tx_info <- get_tx_annotation(range=target_range, tx_data = tx_data)
  
  # Prepare tx info for plotting
  tx_info <- set_tx_level(tx_info, stacking = stacking)
  
  # Make "invisible" transcript plot
  tx_track <- browserly_tx_track(tx_info = tx_info, track_name = 'Annotation')
  
  # Make supblots
  # Decide the type of datatrack to plot if not provided
  if(is.null(type)){
    type <- ifelse(ncol(mcols(cvg)) > hm_thresh, 'heatmap', 'scatter')
  }
  # Convert the GRanges object to a GRanges List
  if(type == 'heatmap') {
    cvg_L = GRangesList(cvg)
    names(cvg_L) <- 'Heatmap'
  } else {
    cvg_L <- sapply(names(mcols(cvg)), function(x) {cvg[, x]})
    }
  
  # Make the subplots
  cvg_tracks <- make_subplots(grl = cvg_L, 
                              type = type, 
                              legend = 'none',
                              ...)
  plots <- c(tx_track, cvg_tracks)
  
  # Set the track heights and plot them
  track_heights <- c(0.2, rep(0.8/length(cvg_tracks), length(cvg_tracks)))
  sp <- subplot(plots, nrows=length(plots), shareX = TRUE, heights=track_heights)
  sp_info <- get_subplot_ax_info(sp)
  
  # Get the axis on which annotations are made
  grt_ax <- get_annotation_axis(ax_info = sp_info)
  
  # Add the annotation shapes
  sp <- add_tx_shapes(plotly_obj = sp, 
                      tx_info = tx_info,
                      target_range = target_range, 
                      grt_ax = grt_ax)
  
  # Adjust the layout parameters
  sp <- modify_y(sp, sp_info, grt_ax, type, sync_y = sync_y)
  
  # Final layout touch up
  sp$x$layout$margin <- list(l = sp$x$layout$margin$l,
                             b = ifelse(is.null(xlabel), 20, 30),
                             t = 10,
                             r = 25)
  sp$x$layout$xaxis$title <- xlabel
  # Force the x axis range to match
  sp$x$layout$xaxis$range <- c(start(target_range), end(target_range))
  return(sp)
  
}

#' Title
#'
#' @param plotly_obj 
#' @param tx_info 
#' @param target_range 
#' @param grt_ax 
#'
#' @return
#' @export
#'
#' @examples
add_tx_shapes <- function(plotly_obj, 
                          tx_info,
                          grt_ax,
                          target_range,
                          y_scaling=0){
  # Add annotations
  tx_info <- biovizBase::mold(tx_info)
  cds_rect <- make_rect(tx_info[tx_info$feature == 'cds', ], height = 0.4, grt_ax)
  utr_rect <- make_rect(tx_info[grep("utr", tx_info$feature), ], height=0.25, grt_ax)
  ncRNA_rect <- make_rect(tx_info[tx_info$feature == 'ncRNA', ], height=0.25, grt_ax)
  
  # Crop arrows to view range
  cropped_introns <- crop_introns(tx_info[tx_info$feature == 'intron', ], target_range)
  intron_arrow <- make_arrows2(cropped_introns, grt_ax, 
                               arrowlen = width(target_range) * 0.01)
  # Compile new shapes
  tx_shapes <- c(cds_rect, utr_rect, ncRNA_rect, intron_arrow)
  
  # Add shapes to list if it already exists
  if(!is.null(plotly_obj$x$layout$shapes)){
    tx_shapes <- c(plotly_obj$x$layout$shapes, tx_shapes)
  }
  plotly_obj$x$layout$shapes <- tx_shapes
  return(plotly_obj)
}



#' browserly_cvg_track
#' Make an interactive coverage track plot that can stand alone or be used as a subplot
#'
#' @param cvg GRanges: mcols should include the data to be plotted
#' @param track_name str: unique identifier for the plot. This is used to find appropriate axes when controlling annotations and domains
#' @param type str: type of plot
#' @param ... additional arguments passed to plotly function "add_trace"
#'
#' @return plotly object
#' @export
#'
#' @examples
browserly_cvg_track <- function(cvg,
                                track_name,
                                type = c('scatter', 'heatmap'),
                                fill = 'tozeroy',
                                mode = 'lines',
                                showlegend = c(TRUE, FALSE),
                                colors = NULL,
                                ...){
  type <- match.arg(type)
  
  # Pull the coverage data as a matrix
  track_data <- t(as.matrix(mcols(cvg)))
  colnames(track_data) <- get_midpoint(cvg)
  
  # Plots are initialized before adding traces. This allows there to be one
  # "reference" plot, which can be used for adding titles or checking domains
  x_data <- as.numeric(colnames(track_data))
  
  if(is.null(colors)){
    colors = RColorBrewer::brewer.pal(6, "Dark2")
    colors = rep(colors, length.out = length(rownames(track_data)))
  }
  
  if(type == 'scatter'){
    # Initialize the plot object that will contain the traces
    p <- plot_ly(source = track_name,
                 type = 'scatter',
                 name = track_name,
                 showlegend = showlegend,
                 colors = colors)
    for(name in rownames(track_data)){
      p <- p %>% add_trace(x = x_data, 
                           y = track_data[name, ],
                           name = name,
                           type = type,
                           hoverinfo='x+y+name',
                           fill = fill,
                           mode = mode,
                           color = name,
                           ...)
    }
  }
  else  if(type == 'heatmap'){
    # A heatmap is only a single trace
    p <- plot_ly(source = track_name,
                 type = 'heatmap',
                 z = track_data,
                 y = rownames(track_data),
                 x=x_data,
                 hoverinfo='x+y+z',
                 name = track_name,
                 colorscale = continuous_colorscale("Purples")(track_data),
                 ...)
  }
  return(p)
}