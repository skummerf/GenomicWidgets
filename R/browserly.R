#' make_browserly_function
#' Abstract browserly plotting to a function that takes a single argument. This is
#' particularly useful when using the function in a shiny app
#'
#' @param cvg_files 
#' @param tx_data 
#' @param type 
#' @param cvg_scaling 
#' @param sample_names 
#' @param stacking 
#' @param sync_y 
#' @param xlabel 
#' @param snps 
#' @param ... 
#'
#' @return
#' @export
#'
#' @author Justin Finkle
#' @examples
make_browserly_function <- function(cvg_files, 
                                    tx_data,
                                    sample_names = names(cvg_files),
                                    cvg_scaling = rep(1, length(cvg_files)),
                                    hm_thresh = 4,
                                    type = NULL,
                                    stacking = c('dense', 'squish'),
                                    sync_y = TRUE,
                                    xlabel = NULL,
                                    snps = GRanges(),
                                    ...){
  plot_browserly <- function(target_range){
    snps_in_range <- get_snps_in_range(snps, target_range)
    print('Getting coverage...')
    cvg_gr <- make_coverage_tracks(inputs = cvg_files,
                                      target_range = target_range, 
                                      sample_names = sample_names, 
                                      scaling_factors = cvg_scaling)
    b_plot <- plot_single_locus(target_range = target_range, 
                                tx_data = tx_data,
                                cvg = cvg_gr,
                                hm_thresh = hm_thresh,
                                type = type,
                                stacking = stacking,
                                sync_y = sync_y,
                                xlabel = xlabel,
                                snps = snps_in_range,
                                ...)
    return(b_plot)
  }
  class(plot_browserly) <- c(class(plot_browserly), "browserly")
  return(plot_browserly)
}


#' plot_multiple_genes
#' Plot several genes as separate tracks, each with multiple traces
#'
#' @param genes character: vector of gene names to be plotted
#' @param centered_cvg list: of GRanges coverage objects that are lined up at the TSS. See `get_centered_gene_info` for more detailes
#' @param centered_tx list: of GRanges annotation objects that are lined up at the TSS. See `get_centered_gene_info` for more detailes
#' @param hm_thresh integer: threshold above which coverage is plotted as a heatmap
#' @param stacking character: how to plot the annotations.
#' @param type character: the plot type. will override hm_thresh
#' @param fill character: fill for plotly traces
#' @param mode character: drawing mode for trace
#' @param sync_y logical: plots subplots on same y scale if TRUE
#' @param xlabel character: title for x axis
#' @param gene_exprs  matrix: gene expression values to be plotted. rownames correspond to genes. columns are different samples and can have duplicate names. colnames are used for putting values in the correct boxplot
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
                                gene_exprs = matrix(),
                                ...){
  # Make transcript tracks
  tx_tracks <- lapply(genes, function(gene){
    tx_info <- centered_tx[[gene]]
    tx_info <- add_tx_stepping(tx_info, stacking = stacking)
    tx_track <- browserly_annotation_track(tx_info = tx_info, 
                                   track_name = paste0(gene,'_Annotation'))
    return(tx_track[[1]])
  })
  
  # Make subplots
  cvg_tracks <- make_subplots(plot_data = cvg_list, 
                              type = type, 
                              fill = fill, 
                              mode = mode,
                              legend = 'first',
                              ...)
  
  # Order tracks, alternating coverage then annotation
  plot_order <- order(rep(seq_along(genes), 2))
  plots <- c(cvg_tracks, tx_tracks)[plot_order]
  sp <- subplot(plots, shareX = TRUE, nrows=length(plots))
  
  # Add expression data if it exists
  if(!all(is.na(gene_exprs))){
    melted_exprs <- melt(gene_exprs)
    plots <- lapply(genes, function(x){
      cur_exprs <- filter(melted_exprs, Var1 == x)
      showlegend <- if(x == genes[[1]]) TRUE else FALSE
      p <- plot_ly(cur_exprs, y=~value, type = 'box', color=~Var2,
                   boxpoints = "all", jitter = 0.3, pointpos=0,
                   showlegend = showlegend, name = x)
      
      p <- p %>% layout(xaxis = list(showticklabels = FALSE))
      return(p)
    })
    
    # Make expression subplot and combine with existing subplot.
    exp_sp <- subplot(plots, nrows=length(genes_entrez))
    sp <- subplot(sp, exp_sp, nrows=1)
    
    # For some reason the yaxis order changes when subplots are made recursively
    # To compensate, change the sort order
    no_exprs <- FALSE
  } else {
    no_exprs <- TRUE
  }
  
  # Get axes information that is used in modifying layouts
  sp_info <- get_subplot_ax_info(sp)
  
  display_range <- c(min(min(start(centered_cvg))), max(max(end(centered_cvg))))
  
  for(gene in genes){
    # Get the axis on which annotations are made
    tx_info <- tx_list[[gene]]
    tx_info <- add_tx_stepping(tx_info, stacking = 'dense')
    ann_ax <- filter(sp_info, grepl(gene, subplot_name), 
                     is_annotation == TRUE)$yaxis[[1]]
    # Add the annotation shapes
    sp <- add_tx_shapes(plotly_obj = sp, 
                        tx_info = tx_info, 
                        target_range = range(tx_info), 
                        ann_ax = ann_ax)
    
  }
  # Adjust the layout parameters
  heights <- rep(c((0.25/length(genes)), (0.75/length(genes))), length(genes))
  sp <- adjust_y_domains(plotly_obj = sp,
                         ax_info = filter(sp_info, type!='box'),
                         heights = heights,
                         sort_decrease = no_exprs)
  sp <- modify_y(plotly_obj = sp, 
                 ax_info = sp_info, 
                 ann_ax = ann_ax, 
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

#' plot_single_locus
#' Plot coverage at a single region. This is the common trackview used
#'
#' @param target_range GRanges: range to be displayed
#' @param tx_data list: annotation information
#' @param cvg GRanges: coverage data. Metadata columns should contain coverage scores
#' @param hm_thresh integer: threshold above which coverage is plotted as a heatmap
#' @param stacking character: how to plot the annotations.
#' @param type character: the plot type. will override hm_thresh
#' @param sync_y logical: plots subplots on same y scale if TRUE
#' @param xlabel character: title for x axis
#' @param snps GRanges: SNPs to display
#' @param ... additional plotly arguments
#'
#' @return
#' @export
#'
#' @examples
plot_single_locus <- function(target_range, 
                              tx_data, 
                              cvg,
                              hm_thresh = 4,
                              stacking = c('dense', 'squish'),
                              type = NULL,
                              sync_y = TRUE,
                              xlabel = NULL,
                              snps = NULL,
                              ...){
  stacking <- match.arg(stacking)
  
  # Make GeneRegion plot
  tx_info <- get_tx_annotation(range=target_range, tx_data = tx_data)
  
  # Prepare tx info for plotting
  tx_info <- add_tx_stepping(tx_info, stacking = stacking)
  
  # Make "invisible" transcript plot
  tx_track <- browserly_annotation_track(tx_info = tx_info, track_name = 'Annotation')
  
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
  
  # Add snps to annotation track
  if(!is.null(snps) & length(snps) > 0){
    annotation_track <- add_snp_to_annotation_track(tx_track[[1]], snps)
  } else {
    annotation_track <- tx_track
  }
  
  # Make the subplots
  cvg_tracks <- make_subplots(plot_data = cvg_L, 
                              type = type, 
                              legend = 'none',
                              ...)
  plots <- c(annotation_track, cvg_tracks)
  
  # Set the track heights and plot them
  track_heights <- c(0.3, rep(0.7/length(cvg_tracks), length(cvg_tracks)))
  sp <- subplot(plots, nrows=length(plots), shareX = TRUE, heights=track_heights)
  sp_info <- get_subplot_ax_info(sp)
  
  # Get the axis on which annotations are made
  ann_ax <- get_annotation_axis(ax_info = sp_info)
  
  # Add the annotation shapes
  sp <- add_tx_shapes(plotly_obj = sp, 
                      tx_info = tx_info,
                      target_range = target_range, 
                      ann_ax = ann_ax)
  
  # Adjust the layout parameters
  sp <- modify_y(sp, sp_info, ann_ax, type, sync_y = sync_y)
  
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
    p <- plot_ly(type = 'scatter',
                 name = track_name,
                 showlegend = showlegend,
                 colors = colors,
                 mode = mode)
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
    p <- plot_ly(type = 'heatmap',
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