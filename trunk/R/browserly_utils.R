#' Title
#'
#' @param gene 
#' @param cvg_files 
#' @param conditions 
#' @param scaling_factor 
#' @param tx_data 
#' @param symbol_table 
#'
#' @return
#' @export
#'
#' @examples
get_centered_gene_info <- function(gene,
                                   extension,
                                   cvg_files,
                                   conditions,
                                   scaling_factor,
                                   tx_data,
                                   symbol_table){
  # Find the gene in the lookup table
  x <- which(symbol_table$symbol == gene)
  
  # Get the gene TSS, range, cvg, and annotation information
  gene_tss <- ifelse(as.character(symbol_table$strand[x]) == "-",
                     symbol_table$end[x],
                     symbol_table$start[x])
  gene_chr <- as.character(symbol_table$chr[x])
  gene_strand <- as.character(symbol_table$strand[x])
  gene_range <- get_view_range(chr = gene_chr,
                               start = gene_tss - extension,
                               end = gene_tss + extension,
                               strand = gene_strand)
  gene_cvg <- make_coverage_tracks(cvg_files,
                                   target_range = gene_range,
                                   sample_names = conditions,
                                   scaling_factors = scaling_factor)
  gene_tx <- get_tx_annotation(range=gene_range, tx_data = tx_data)
  
  # Center everything around the TSS
  gene_cvg <- shift(gene_cvg, shift = -gene_tss)
  
  # Remove extra seqnames in annotation to allow shifting
  gene_tx <- keepSeqlevels(gene_tx, gene_chr)
  seqlengths(gene_tx) <- NA
  gene_tx <- shift(gene_tx, shift = -gene_tss)
  
  # If the reference gene is on the - strand, the info needs to be flipped
  if(gene_strand == '-'){
    gene_cvg <- invert_ranges(gene_cvg)
    gene_tx <- invert_ranges(gene_tx)
  }
  
  return(list(gene_cvg = gene_cvg, gene_tx = gene_tx))
}

# =============================================================================
# =============================================================================
# Helper functions - not exported
# =============================================================================
# =============================================================================

#' invert_ranges
#' Flip the range over a reference value
#' 
#' @param gr 
#' @param reference 
#'
#' @return
#'
#' @examples
invert_ranges <- function(gr,
                          reference = 0){
  # Convert to a dataframe so GRanges doesn't complain about widths during transformation
  old_start <- start(gr)
  names(gr) <- NULL
  tmp_df <- as.data.frame(gr)
  tmp_df$start <- reference - tmp_df$end
  tmp_df$end <- reference - old_start
  
  return(GRanges(tmp_df))
}

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
                               target_range = NULL,
                               opacity = 0,
                               markercolor = 'grey'){
  # Make GRanges into df
  tx_df <- biovizBase::mold(tx_info) %>% group_by(feature)
  tx_df['text'] <- paste0('Tx ID: ', tx_df$transcript)
  
  # Plot the annotation features
  p <- plot_ly(tx_df, 
               type='scatter',
               name = track_name) %>% 
    add_markers(x = ~midpoint, 
                y=~stepping, 
                showlegend=FALSE,
                text=~text, 
                marker = list(color = markercolor),
                hoverinfo='x+name+text',
                opacity = opacity)
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
                          legend = c('first', 'all', 'none'),
                          ...){
  legend = match.arg(legend)
  
  # Make the subplots
  plot_list <- lapply(names(grl), function(x, grl, type, ...){
    if(x == names(grl)[[1]] & legend == 'first'){
      showlegend = TRUE
    } else if(legend == 'all'){
      showlegend = TRUE
    } else {
      showlegend = FALSE
    }
    browserly_cvg_track(cvg = grl[[x]],
                        track_name = x,
                        type = type,
                        showlegend = showlegend,
                        ...)
  }, grl, type, ...)
  return(plot_list)
}

#' modify_y
#' Modify y axes features such as titles, range, etc
#'
#' @param plotly_obj 
#' @param ax_info 
#' @param grt_ax 
#' @param type 
#' @param sync_y 
#' @param title_rotation 
#'
#' @return
#'
#' @examples
modify_y <- function(plotly_obj, 
                     ax_info, 
                     grt_ax, 
                     type,
                     sync_y = TRUE,
                     native_title = FALSE,
                     title_rotation = 0){
  
  # Get the data traces and their axes
  plot_axes <- filter(ax_info, is_trace == FALSE | type == 'heatmap')
  
  # Get the max y value
  score_max <- max(unlist(filter(ax_info, is_trace == TRUE, yaxis != grt_ax)$ymax))
  
  # Get the largest title and pad with y axis labels
  label_padding <- nchar(round(score_max))+3 # Give a charater buffer in case of zooming
  max_title <- calc_title_margin(plot_axes$subplot_name)
  
  title_list <- list()
  # Adjust the plot parameters for each subplot
  for(idx in 1:nrow(plot_axes)){
    # Get the row as a named list
    cur_row <- lapply(slice(plot_axes, idx), function(x) {unlist(x)})
    ax <- cur_row$yaxis
    if(!cur_row$is_annotation){
      if(cur_row$type!='heatmap'){
        if(sync_y) plotly_obj$x$layout[[ax]][['range']] <- c(0, score_max)
        if(native_title) plotly_obj $x$layout[[ax]][['title']] <- cur_row$subplot_name
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
    
    # Set the title for the axis if necessary
    if(cur_row$type!='heatmap' & !native_title){
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
      # Set the margin and pad for the size of the y axes labels
      plotly_obj$x$layout$margin$l <- max_title
      
      # Add the titles
      plotly_obj$x$layout$annotations <- title_list
    }
  }
  
  return(plotly_obj)
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
                      is_trace = is_trace,
                      is_annotation = grepl('Annotation', x$name))
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
  ann_info <- ax_info$yaxis[grepl(annotation_str, ax_info$subplot_name)]
  if(length(ann_info) == 0){
    stop("No annotation track found in axes information")
  }
  return(unlist(ann_info))
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