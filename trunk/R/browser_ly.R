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
#' @param txdb 
#' @param range 
#' @param tx_data 
#' @param cvg
#' @param hm_thresh 
#' @param type 
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
  cvg_tracks <- make_subplots(grl = cvg_L, type = type, ...)
  plots <- c(tx_track, cvg_tracks)
  
  # Set the track heights and plot them
  track_heights <- c(0.3, rep(0.7/length(cvg_tracks), length(cvg_tracks)))
  sp <- subplot(plots, nrows=length(plots), shareX = TRUE, heights=track_heights)
  sp_info <- get_subplot_ax_info(sp)
  
  # Get the axis on which annotations are made
  grt_ax <- get_annotation_axis(ax_info = sp_info)
  
  # Add the annotation shapes
  sp <- add_tx_shapes(sp, tx_info, target_range, grt_ax)
  
  # Adjust the layout parameters
  sp <- modify_y(sp, sp_info, grt_ax, type, sync_y = sync_y)
  
  # Force the x axis range to match
  sp$x$layout$xaxis$range <- c(start(target_range), end(target_range))
  return(sp)
  
}

modify_y <- function(plotly_obj, 
                     ax_info, 
                     grt_ax, 
                     type,
                     sync_y = TRUE,
                     title_rotation = 0){
  
  # Get the data traces and their axes
  plot_axes <- filter(ax_info, is_trace == FALSE | type == 'heatmap')
  
  # Get the max y value
  score_max <- max(unlist(filter(ax_info, is_trace == TRUE, yaxis != grt_ax)$ymax))
  
  # Get the largest title and pad with y axis labels
  label_padding <- nchar(round(score_max))+1
  max_title <- calc_title_margin(plot_axes$subplot_name)
  
  title_list <- list()
  # Adjust the plot parameters for each subplot
  for(idx in 1:nrow(plot_axes)){
    # Get the row as a named list
    cur_row <- lapply(slice(plot_axes, idx), function(x) {unlist(x)})
    ax <- cur_row$yaxis
    if(ax != grt_ax){
      if(cur_row$type!='heatmap'){
        if(sync_y) plotly_obj$x$layout[[ax]][['range']] <- c(0, score_max)
      } else {
        # Get tick labels to update margin size
        hm_idx <- unlist(filter(ax_info, yaxis == ax, is_trace == TRUE)$data_idx)
        # Use the index of the current heatmap to pull the y values from the plotly data structure
        # which correspond to the tick names
        max_title <- max(max_title, calc_title_margin(plotly_obj$x$data[[hm_idx]]$y))
        plotly_obj$x$layout[[ax]]['showline'] <- TRUE
        plotly_obj$x$layout[[ax]]['mirror'] <- TRUE
        plotly_obj$x$layout[[ax]]['ticks'] <- ""
      }
    } else {
      # Modify transcript axis
      plotly_obj$x$layout[[ax]][['ticks']] <- ""
      plotly_obj$x$layout[[ax]][['showticklabels']] <- FALSE
      plotly_obj$x$layout[[ax]][['showgrid']] <- FALSE
    }
    
    # Set the title for the axis
    if(cur_row$type!='heatmap'){
      title_list[[idx]] <- list(text = cur_row$subplot_name, 
                                showarrow=FALSE, 
                                textangle=title_rotation,
                                x = cur_row$x0-(label_padding/100),
                                y = (cur_row$y0+cur_row$y1)/2,
                                xref = 'paper',
                                yref = 'paper',
                                borderpad = 0,
                                xanchor = 'right'
      )
    }
  }
  # Set the margin and pad for the size of the y axes labels
  plotly_obj$x$layout$margin$l <- max_title
  
  # Add the titles
  plotly_obj$x$layout$annotations <- title_list
  
  return(plotly_obj)
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
add_tx_shapes <- function(plotly_obj, tx_info, target_range, grt_ax, y_scaling=0){
  # Add annotations
  tx_info <- biovizBase::mold(tx_info)
  cds_rect <- make_rect(tx_info[tx_info$feature == 'cds', ], height = 0.4, grt_ax)
  utr_rect <- make_rect(tx_info[grep("utr", tx_info$feature), ], height=0.25, grt_ax)
  ncRNA_rect <- make_rect(tx_info[tx_info$feature == 'ncRNA', ], height=0.25, grt_ax)
  
  # Crop arrows to view range
  cropped_introns <- crop_introns(tx_info[tx_info$feature == 'intron', ], target_range)
  intron_arrow <- make_arrows2(cropped_introns, grt_ax, 
                               arrowlen = width(target_range) * 0.01)
  tx_shapes <- c(cds_rect, utr_rect, ncRNA_rect, intron_arrow)
  tx_shapes <- lapply(tx_shapes, function(x, y_scaling){
    if(x$type == 'rect'){
      x$y0 <- x$y0-y_scaling*0.1
    } else if(x$type == 'line'){
      x$y0 <- x$y0-y_scaling*0.05
      x$y1 <- x$y1-y_scaling*0.05
    }
    
    return(x)
  }, y_scaling)
  if(!is.null(plotly_obj$x$layout$shapes)){
    tx_shapes <- c(plotly_obj$x$layout$shapes, tx_shapes)
  }
  grt_layout <- list(showlegend=FALSE)
                     #annotations=intron_arrow)
  plotly_obj$x$layout <- modifyList(plotly_obj$x$layout, grt_layout)
  plotly_obj$x$layout$shapes <- tx_shapes
  # grt_y_layout <- list(autorange='reversed', showticklabels=FALSE, showticks=FALSE,
  #                      title='Transcripts')
  # plotly_obj$x$layout[[grt_ax]] <- modifyList(plotly_obj$x$layout[[grt_ax]], grt_y_layout)
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
                                ...){
  type <- match.arg(type)
  
  # Pull the coverage data as a matrix
  track_data <- t(as.matrix(mcols(cvg)))
  colnames(track_data) <- get_midpoint(cvg)
  
  # Plots are initialized before adding traces. This allows there to be one
  # "reference" plot, which can be used for adding titles or checking domains
  x_data <- as.numeric(colnames(track_data))
  if(type == 'scatter'){
    # Initialize the plot object that will contain the traces
    p <- plot_ly(source = track_name,
                 type = 'scatter',
                 name = track_name)
    for(name in rownames(track_data)){
      p <- p %>% add_trace(x = x_data, 
                           y = track_data[name, ],
                           name = name,
                           type = type,
                           hoverinfo='x+y+name',
                           fill = fill,
                           mode = mode,
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



# =============================================================================
# =============================================================================
# Helper functions - not exported
# =============================================================================
# =============================================================================

#' browserly_tx_track
#' Make the invisible points for an annotation track.
#' 
#' Currently this is only a helper function, but it could be made to stand alone
#'
#' @param tx_info GRanges: mcols describe the annotation features
#'
#' @return plotly_object
#'
#' @examples
browserly_tx_track <- function(tx_info, 
                               track_name,
                               add_shapes = FALSE,
                               target_range = NULL){
  # Make GRanges into df
  tx_df <- biovizBase::mold(tx_info) %>% group_by(feature)
  tx_df['text'] <- paste0('Tx ID: ', tx_df$transcript)
  
  # Plot the annotation features
  p <- plot_ly(tx_df, 
               type='scatter',
               name = track_name) %>% 
    add_markers(x = ~midpoint, 
                y=~stepping, 
                color=~feature, 
                showlegend=FALSE,
                text=~text, 
                hoverinfo='x+name+text',
                opacity = 0)
  target_range <- if(is.null(target_range)) reduce(range(tx_info), 
                                                   ignore.strand=TRUE)
  # Build the plot so it can be referenced. Return as list for subploting
  p <- list(plotly_build(p))
  return(p)
}

#' make_subplots
#' Make subplots of coverage tracks
#'
#' @param grl GRangesList: a named list. Each GRanges object in the list will be plotted on a separate axis
#' @param type str: type of plot to use
#' @param ... additional parameters to pass to browserly_cvg_track
#'
#' @return
#'
#' @examples
make_subplots <- function(grl,
                          type,
                          ...){
  # Make the subplots
  plot_list <- lapply(names(grl), function(x, grl, type, ...){
    browserly_cvg_track(cvg = grl[[x]],
                        track_name = x,
                        type = type,
                        ...)
  }, grl, type, ...)
  return(plot_list)
}

#' get_subplot_ax_info
#' Get the axes information in a plotly object. This is meant to be
#' called on an object created by the subplot function
#'
#' @param plotly_obj 
#'
#' @return
#'
#' @examples
get_subplot_ax_info <- function(plotly_obj){
  sp_ax_info <- sapply(seq_along(plotly_obj$x$data), function(idx){
    # Get the plot info
    # If there is x, y, or z data then it is a trace
    x <- plotly_obj$x$data[[idx]]
    is_trace <- (is.numeric(x$x) | is.numeric(x$y) | is.numeric(x$z))
    data_info <- list(subplot_name = x$name,
                      data_idx = idx,
                      xid = x$xaxis,
                      yid = x$yaxis,
                      xaxis = gsub("x", "xaxis", x$xaxis),
                      yaxis = gsub("y", "yaxis", x$yaxis),
                      ymax = ifelse(is.numeric(x$y), max(x$y), 0),
                      type = x$type,
                      is_trace = is_trace)
    return(data_info)
  })
  
  # Get the domains for the axes
  domain_info <- apply(sp_ax_info, 2, function(x, plotly_obj) {
    # Comine the domains
    cbind(plotly_obj$x$layout[[x$yaxis]]$domain,
          plotly_obj$x$layout[[x$xaxis]]$domain)
    }, plotly_obj)
  rownames(domain_info) <- c('y0', 'y1', 'x0', 'x1')
  sp_ax_info <- rbind(sp_ax_info, domain_info)
  sp_ax_info <- data.frame(t(sp_ax_info))
  return(sp_ax_info)
}


#' calc_title_margin
#' Calculate how many pixels of space the titles need. Assumes the titles are horizontal
#'
#' @param titles vector
#' @param letter_width int: number of pixels each letter requires
#' @param padding int: number of character to pad to the maximium title found
#'
#' @return
#'
#' @examples
calc_title_margin <- function(titles,
                              padding = 0,
                              letter_width = 9){
  # Calculate how big the left margin should be
  max_label_length <- max(sapply(titles, nchar)) + padding
  left_margin = max(letter_width*max_label_length, 60)
  return(left_margin)
}


#' get_annotation_axis
#'
#' @param ax_info data.frame like: contains info for different subplot axes in plotly object
#' @param annotation_str str: the string used to find the annotation axis
#'
#' @return str: the annotation axis reference in the form "yaxis#"
#'
#' @examples
get_annotation_axis <- function(ax_info,
                                annotation_str = 'Annotation'){
  ann_info <- filter(ax_info, subplot_name == annotation_str)$yaxis
  if(length(ann_info) == 0){
    stop("No annotation track found in axes information")
  }
  ann_ax <- ann_info[[1]]
  if(length(ann_info) > 1){
    warning(paste0('Multiple possible annotation axes were found.',
                   'Using the first one encountered: ', ann_ax))
  }
  return(ann_ax)
}

#' set_tx_level
#' Set the level that that transcripts will be displayed on
#'
#' @param tx_gr GRanges: contains transcript information 
#' @param stacking str: dense collapses transcripts, squish expands them
#'
#' @return
#'
#' @examples
set_tx_level <- function(tx_gr, 
                         stacking = c('dense', 'squish')){
  stacking <- match.arg(stacking)
  if(length(tx_gr)){
    if(stacking == 'squish'){
      tx_gr <- biovizBase::addStepping(tx_gr, group.name = "transcript",
                                       group.selfish = FALSE)
    } else if(stacking == 'dense'){
      # The simple solution for now
      mcols(tx_gr)$stepping <- 1
    }
    return(tx_gr)
  }
}

crop_introns <- function(introns, target_range){
  introns$start <- pmax(introns$start, start(target_range))
  introns$end <- pmin(introns$end, end(target_range))
  introns$midpoint <- pmin(pmax(introns$midpoint, start(target_range)), end(target_range))
  return(introns)
}

make_rect <- function(df, height, yref){
  yref <- gsub("yaxis", "y", yref)
  if(nrow(df)>0){
    rect_list <- vector("list", nrow(df))
    for(e in 1:nrow(df)){
      row <- df[e, ]
      rect_list[[e]] <- list(type = "rect", fillcolor = "blue", opacity = 1, line=list(width=0),
                             x0 = row$start, x1 = row$end, xref = "x",
                             y0 = row$stepping-height, y1 = row$stepping+height, yref = yref)
    }
  } else {
    rect_list <- NULL
  }
  return(rect_list)
}

make_arrows <- function(df, yref){
  yref <- gsub("yaxis", "y", yref)
  if(nrow(df)>0){
    arrow_list <- vector("list", nrow(df))
    for(i in 1:nrow(df)){
      row <- df[i, ]
      xstart <- ifelse(row$strand =="-", row$start, row$end)
      xend <- ifelse(row$strand =="-", row$end, row$start)
      arrow_list[[i]] <- list(x = xstart, y=row$stepping, xref = "x", yref = yref,
                              showarrow = TRUE, ax = xend, ay=row$stepping,
                              axref='x', ayref= yref, arrowwidth = 1, text="")
      
    }
  } else {
    arrow_list <- NULL
  }
  return(arrow_list)
}

arrow_helper <- function(arrow_start, strand, arrowlen, arrowheight, y, yref){
  arrow_end <- ifelse(strand =="-", arrow_start + arrowlen, arrow_start - arrowlen)
  list(list(x0 = arrow_start, x1=arrow_end, 
            y0 = y, y1 = y - arrowheight, 
            xref = "x", yref = yref,
            type = "line",
            line = list(width = 0.5)),
       list(x0 = arrow_start, x1=arrow_end, 
            y0 = y, y1 = y + arrowheight, 
            xref = "x", yref = yref,
            type = "line",
            line = list(width = 0.5)))
}

make_arrows2 <- function(df, yref, arrowlen = 500, arrowheight = 0.15, arrowgap = 1500){
  yref <- gsub("yaxis", "y", yref)
  if(nrow(df)>0){
    line_list <- vector("list", nrow(df))
    arrow_list <- vector("list", nrow(df))
    for(i in 1:nrow(df)){
      row <- df[i, ]
      xstart <- 
        xend <- ifelse(row$strand =="-", row$end, row$start)
      line_list[[i]] <- list(x0 = row$start, y0=row$stepping, 
                             x1 = row$end, y1 = row$stepping, 
                             xref = "x", yref = yref,
                             type = "line",
                             line = list(width = 0.5, 
                                         dash = ifelse(row$strand == "-", "dot","solid")))
      if (row$end - row$start > 2 * arrowlen){
        if (row$strand == "-"){
          arrow_pos <- row$midpoint - arrowlen * 0.5
        } else{
          arrow_pos <- row$midpoint + arrowlen * 0.5
        }
        arrow_list[[i]] <- arrow_helper(arrow_pos,
                                        row$strand,
                                        arrowlen,
                                        arrowheight,
                                        row$stepping,
                                        yref)
      }
    }
    out <- c(unlist(arrow_list, recursive = FALSE), line_list)
  } else {
    out <- NULL
  }
  return(out)
}
# =============================================================================
# =============================================================================
# Deprecated Functions
# =============================================================================
# =============================================================================

coverage_heatmap <- function(cvg){
  hm_vals <- apply(t(as.matrix((mcols(cvg)))), 2, rev)
  midpoint <- get_midpoint(cvg)
  p <- plot_ly(z=hm_vals, y=rev(colnames(mcols(cvg))), x=midpoint, type='heatmap',
               hoverinfo='x+z', colorscale = continuous_colorscale("Purples")(hm_vals))
  return(p)
}

collapse_tx <- function(gr){
  # Reduce introns and RNA into minimal set of ranges
  introns <- gr[mcols(gr)$feature == 'intron']
  dense_introns <- reduce(introns)
  mcols(dense_introns)$feature <- 'intron'
  rna <- gr[mcols(gr)$feature != 'intron']
  dense_rna <- reduce(rna)
  
  # Add metadata to RNA to keep track of transcripts and features
  overlaps <- findOverlaps(rna, dense_rna)
  mcols(dense_rna)$feature <- splitAsList(mcols(rna)$feature[queryHits(overlaps)],
                                          factor(subjectHits(overlaps)))
  mcols(dense_rna)$transcript <- splitAsList(mcols(rna)$transcript[queryHits(overlaps)],
                                             factor(subjectHits(overlaps)))
  return(c(dense_introns, dense_rna))
}
