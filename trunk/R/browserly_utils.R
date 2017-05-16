# =============================================================================
# =============================================================================
# Exported functions
# =============================================================================
# =============================================================================


#' get_centered_gene_info
#' Center coverage and annotation around TSS
#' NOTES: 1. Can be generalized to pass target ranges instead of gene list
#'        2. Can be generalized to center around any feature
#'        3. Probably should only pass annotations relevant to the genes
#'        
#' @param gene character: gene symbol to be centered
#' @param extension integer: upstream (and downstream) region around center to 
#' be viewed
#' @param cvg_files characters: paths to coverage files of interest
#' @param sample_names characters: names for the cvg_files to be displayed
#' @param scaling_factor numeric: scaling factor for coverage, e.g. library size
#' @param tx_data list: annotation information
#' @param symbol_table data.frame: annotation information linking gene symbols 
#' with txdb information
#'
#' @return
#' @export
#' @author Justin Finkle
#'
#' @examples
get_centered_gene_info <- function(gene,
                                   extension = 50000,
                                   cvg_files,
                                   sample_names,
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
  gene_range <- extend_grange(gr = get_view_range(chr = gene_chr,
                               start = gene_tss,
                               end = gene_tss,
                               strand = gene_strand),
                              extend = extension)
  gene_cvg <- make_coverage_tracks(cvg_files,
                                   target_range = gene_range,
                                   sample_names = sample_names,
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
    gene_cvg <- reflect_ranges(gene_cvg)
    gene_tx <- reflect_ranges(gene_tx)
    strand(gene_tx) <- "*"
  }
  
  return(list(gene_cvg = gene_cvg, gene_tx = gene_tx))
}

#' get_snps_in_range
#' Get the SNPs that fall in a particluar range
#'
#' @param snp_gr GRanges: contains SNP information, e.g. from GWAS catalog. See 
#' load_snp_data for more info
#' @param target_range GRanges: the range in which to look for SNPs
#'
#' @return
#' @export
#'
#' @author Justin Finkle
#' @examples
get_snps_in_range <- function(snp_gr = GRanges(), 
                              target_range = GRanges()){
  
  snp_hits <- findOverlaps(snp_gr, target_range, type='within')
  snp_locs <- snp_gr[queryHits(snp_hits)]
  
  return(snp_locs)
}

# =============================================================================
# =============================================================================
# Helper functions - not exported
# =============================================================================
# =============================================================================

#' adjust_y_domains
#' Meant to adjust yaxes domains, primarily for mesocale views. Adjusts domains
#' to properly view annotations with each subplot
#'
#' @param plotly_obj plotly: built plotly object to be modified
#' @param ax_info data.frame: contains relevant axes information. created using 
#' get_subplot_ax_info
#' @param heights numeric: heights of each track, including annotation tracks
#'
#' @return plotly_obj
#'
#' @author Justin Finkle
#' @examples
adjust_y_domains <- function(plotly_obj, 
                             ax_info, 
                             heights,
                             sort_decrease = TRUE){
  # The sort order may change depending on external factors, such as expression 
  # data
  yaxes <- sort(unlist(unique(ax_info$yaxis)), decreasing = sort_decrease)
  tops <- cumsum(heights)
  bases <- tops - heights
  
  # Check that heights sum to less than 1 and there are the same number of axes 
  # and heights
  stopifnot(sum(heights) <= 1, length(yaxes) == length(heights))
  
  for(y_idx in seq_along(yaxes)){
    yax <- yaxes[[y_idx]]
    plotly_obj$x$layout[[yax]]$domain <- c(bases[[y_idx]], tops[[y_idx]])
  }
  return(plotly_obj)
}

#' reflect_ranges
#' Flip the range over a reference value
#' 
#' @param gr GRanges: data to be flipped
#' @param reference integer: position around which values in the gr are reflected
#'
#' @return
#'
#' @author Justin Finkle
#' @examples
reflect_ranges <- function(gr,
                          reference = 0){
  # Convert to a dataframe so GRanges doesn't complain about widths during 
  # transformation
  old_start <- start(gr)
  
  # Remove names so there isn't a conflict when coercing to a data.frame
  names(gr) <- NULL
  tmp_df <- as.data.frame(gr)
  tmp_df$start <- reference - tmp_df$end
  tmp_df$end <- reference - old_start
  
  return(GRanges(tmp_df))
}

#' add_snp_to_annotation_track
#' Add SNP information to an annotation track.
#' NOTE: SNPS could be made as a separate track. However, additional subplots 
#' compress
#' the vertical space, and when there is insufficient room the hoverinfo goes 
#' away.
#'
#' @param ann_track plotly: annotation track
#' @param snp_info GRanges: contains SNP information. Currently supports GRanges
#'  made from GWAS catalog
#' @param group_col character: the column in snp_info for grouping traces. 
#' Default is "CONTEXT" which corresponds to the type of mutation in the GWAS 
#' catalog, e.g. 'missense'
#'
#' @return plotly_object
#'
#' @author Justin Finkle
#' @examples
add_snp_to_annotation_track <- function(ann_track,
                                        snp_info, 
                                        group_col = 'CONTEXT'){
  # Make GRanges into df and group into traces.
  snp_df <- biovizBase::mold(snp_info) %>% group_by(.dots = c(group_col))
  
  # Make the hoverinfo text
  snp_df['text'] <- paste0('SNP ID: ', names(snp_info), "<br>",
                           'Disease: ', snp_df$DISEASE.TRAIT, "<br>",
                           'Risk Allele (freq): ', 
                           sapply(snp_df$STRONGEST.SNP.RISK.ALLELE, 
                                  function(x){substr(x, nchar(x), nchar(x))}), 
                           " (",snp_df$RISK.ALLELE.FREQUENCY, ")")
  
  # Set the level of the snps. Each type of SNP gets a different stepping.
  # These are forced to be negative so they appear below annotations
  snp_types <- unique(snp_df[[group_col]])
  snp_df['level'] <- vapply(snp_df[[group_col]], 
                            function(x){-which(x==snp_types)}, 
                            FUN.VALUE = 0L)

  p <- ann_track %>% add_markers(x = snp_df[['midpoint']],
                                y = snp_df[['level']],
                                color = snp_df[[group_col]],
                                showlegend=FALSE,
                                text=snp_df[['text']],
                                hoverinfo='x+text')
  
  # Build the plot so it can be referenced. Return as list for subploting
  p <- list(plotly_build(p))
  return(p)
}



#' browserly_annotation_track
#' Make the invisible points for an annotation track.
#' 
#' Currently this is only a helper function, but it could be made to stand alone
#'
#' @param tx_info GRanges: mcols describe the annotation features
#' @param track_name character: name for the track. May be used for referencing 
#' axes and adjusting domains
#' @param opacity numeric: value between 0 and 1. Sets the opacity of the points 
#' used for hoverinfo. Should almost always be 0.
#' @param markercolor character: color to use for markers. Default is 'grey'. If 
#' set to NULL, different colors will be used for each class of annotation 
#' feature, e.g. intron, cds, utr3
#'
#' @return plotly_object
#'
#' @author Justin Finkle
#' @examples
browserly_annotation_track <- function(tx_info, 
                                       track_name = 'Annotation',
                                       opacity = 0,
                                       markercolor = 'grey'){
  # Convert the marker color to a plotly usable color
  markercolor <- make_plotly_color(markercolor)
  
  if (length(tx_info) > 0){
    # Make GRanges into df and group by annotation feature
    tx_df <- biovizBase::mold(tx_info) %>% group_by(feature)
    tx_df['text'] <- paste0('Tx ID: ', tx_df$transcript)
   
    # Plot the annotation features
    p <- plot_ly(tx_df, 
               type='scatter',
               name = track_name,
               mode = 'markers') %>% 
      add_markers(x = ~midpoint, 
                y=~stepping,
                color = ~feature,
                colors = brewer.pal.helper(length(unique(tx_df$feature)), 
                                           "Set2"),
                showlegend = FALSE,
                text = ~text, 
                marker = list(color = markercolor),
                hoverinfo = 'x+name+text',
                opacity = opacity)
  } else{
    p <- plot_ly(type='scatter',
                 name = track_name,
                 mode = 'markers')
  }
  
  # Build the plot so it can be referenced. Return as list for subploting
  p <- list(plotly_build(p))
  return(p)
}

#' make_subplots
#' Make subplots of coverage tracks
#' 
#'
#' @param plot_data list: a named list. Each list item will be plotted on a 
#' separate axis, and is expected to be a GRanges object, whose metadata columns
#'  will be plotted.
#' @param type str: type of plot to use
#' @param legend str: which subplots to include with legends. 'first' (default), 
#' only shows the legend for the first subplot.
#' @param ... additional parameters to pass to browserly_cvg_track
#'
#' @return
#'
#' @author Justin Finkle
#' @examples
make_subplots <- function(plot_data,
                          type,
                          legend = c('first', 'all', 'none'),
                          ...){
  legend = match.arg(legend)
  
  # Make the subplots
  plot_list <- lapply(names(plot_data), function(x, plot_data, type, ...){
    if(x == names(plot_data)[[1]] & legend == 'first'){
      showlegend = TRUE
    } else if(legend == 'all'){
      showlegend = TRUE
    } else {
      showlegend = FALSE
    }
    browserly_cvg_track(cvg = plot_data[[x]],
                        track_name = x,
                        type = type,
                        showlegend = showlegend,
                        ...)
  }, plot_data, type, ...)
  return(plot_list)
}

#' modify_y
#' Modify y axes features such as titles, range, etc. This is a HIGHLY specific
#' function for making trackviews look nice.
#'
#' @param plotly_obj plotly: object to modify
#' @param ax_info data.frame: contains axes information so plot can be adjusted. 
#' See get_subplot_ax_info for more information
#' @param ann_ax character: the annotation axis of the form "yaxis[2-9]". 
#' These have different formatting requirements
#' @param type character: different types of plots are treated differently
#' @param sync_y logical: put all y-axes on same scale
#' @param native_title logical: use plotly titles. If FALSE, annotations are 
#' added to each subplot. Useful if titles are long to prevent title overlap.
#'  Margins are automatically adjusted
#' @param title_rotation numeric: angle of plot titles. Only used if 
#' native_title is FALSE
#'
#' @return
#'
#' @author Justin Finkle
#' @examples
modify_y <- function(plotly_obj, 
                     ax_info, 
                     ann_ax, 
                     type,
                     sync_y = TRUE,
                     native_title = FALSE,
                     title_rotation = 0){
  
  # Get the data traces and their axes
  plot_axes <- filter(ax_info, is_trace == FALSE | type == 'heatmap')
  
  # Get the max y value
  score_max <- max(unlist(filter(ax_info, is_trace == TRUE, 
                                 yaxis != ann_ax)$ymax))
  
  # Get the largest title and pad with y axis labels
  # Give a charater buffer in case of zooming
  label_padding <- nchar(round(score_max))+3 
  max_title <- calc_title_margin(plot_axes$subplot_name)
  
  title_list <- list()
  # Adjust the plot parameters for each subplot
  for(idx in 1:nrow(plot_axes)){
    # Get the row as a named list
    cur_row <- lapply(dplyr::slice(plot_axes, idx), function(x) {unlist(x)})
    ax <- cur_row$yaxis
    if(!cur_row$is_annotation & !cur_row$is_snp){
      if(cur_row$type!='heatmap'){
        if(sync_y) plotly_obj$x$layout[[ax]][['range']] <- c(0, score_max)
        if(native_title) plotly_obj$x$layout[[ax]][['title']] <- 
            cur_row$subplot_name
      } else {
        # Get tick labels to update margin size
        hm_idx <- unlist(filter(ax_info, yaxis == ax, 
                                is_trace == TRUE)$data_idx)
        # Use the index of the current heatmap to pull the y values from the 
        # plotly data structure which correspond to the tick names
        max_title <- max(max_title, 
                         calc_title_margin(plotly_obj$x$data[[hm_idx]]$y))
        plotly_obj$x$layout[[ax]]['ticks'] <- ""
        
        # Add a box around the heatmap
        box <- list(list(type = 'rect', x0 = cur_row$x0, x1 = cur_row$x1,
                    y0 = cur_row$y0, y1 = cur_row$y1,
                    xref= 'paper', yref='paper', line = list(width=1)))
        plotly_obj$x$layout$shapes <- append(plotly_obj$x$layout$shapes, box)
      }
    } else {
      # Modify transcript axis
      plotly_obj$x$layout[[ax]][['ticks']] <- ""
      plotly_obj$x$layout[[ax]][['showticklabels']] <- FALSE
      plotly_obj$x$layout[[ax]][['showgrid']] <- FALSE
      plotly_obj$x$layout[[ax]][['zeroline']] <- FALSE
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
      )}
  }
  # Set the margin and pad for the size of the y axes labels
  plotly_obj$x$layout$margin$l <- max_title + 2*label_padding
  
  # Add the titles
  plotly_obj$x$layout$annotations <- title_list
  
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
#' @author Justin Finkle
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
                      is_annotation = grepl('Annotation', x$name),
                      is_snp = grepl('SNP', x$name))
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
#' Calculate how many pixels of space the titles need. Assumes the titles are 
#' horizontal
#'
#' @param titles vector
#' @param letter_width int: number of pixels each letter requires
#' @param padding int: number of character to pad to the maximium title found
#'
#' @return
#'
#' @author Justin Finkle
#' @examples
calc_title_margin <- function(titles,
                              padding = 0,
                              letter_width = 9){
  # Calculate how big the left margin should be
  max_label_length <- max(vapply(titles, nchar, FUN.VALUE = 0L)) + padding
  left_margin = max(letter_width*max_label_length, 60)
  return(left_margin)
}


#' get_annotation_axis
#'
#' @param ax_info data.frame like: contains info for different subplot axes in 
#' plotly object
#' @param annotation_str str: the string used to find the annotation axis
#'
#' @return str: the annotation axis reference in the form "yaxis[2-9]"
#'
#' @author Justin Finkle
#' @examples
get_annotation_axis <- function(ax_info,
                                annotation_str = 'Annotation'){
  ann_info <- ax_info$yaxis[grepl(annotation_str, ax_info$subplot_name)]
  if(length(ann_info) == 0){
    stop("No annotation track found in axes information")
  }
  return(unlist(ann_info))
}

#' add_tx_stepping
#' Set the level that that transcripts will be displayed on
#'
#' @param tx_gr GRanges: contains transcript information 
#' @param stacking str: dense collapses transcripts, squish expands them
#'
#' @return
#'
#' @author Justin Finkle
#' @examples
add_tx_stepping <- function(tx_gr, 
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

#' make_plotly_color
#' Convert R colors to plotly compatible rgb values
#'
#' @param color_str 
#'
#' @return
#' @export
#'
#' @author Justin Finkle
#' @examples
make_plotly_color <- function(color_str){
  # Plotly accepts some string colors, but it is safer to use rgb(0,0,0) values
  p_color <- rgb(t(col2rgb(color_str)), maxColorValue = 255)
  return(p_color)
}

# =============================================================================
# =============================================================================
# Annotation track shapes
# =============================================================================
# =============================================================================
#' add_tx_shapes
#' Add annotation shapes to a specific axis
#'
#'
#' @param plotly_obj plotly: object to be modified
#' @param tx_info GRanges: contains annotation information. Must have metadata column with name 'feature' to make shapes
#' @param ann_ax character: axis to add shapes to, of the form 'yaxis[2-9]'
#' @param target_range GRanges: range to display. Used to draw shapes properly
#'
#' @return
#' @export
#'
#' @author Justin Finkle
#' @examples
add_tx_shapes <- function(plotly_obj, 
                          tx_info,
                          ann_ax,
                          target_range){
  
  if (length(tx_info) > 0){
    
    # Add annotations
    tx_info <- biovizBase::mold(tx_info)
    cds_rect <- make_rect(tx_info[tx_info$feature == 'cds', ],
                          height = 0.4, 
                          ann_ax)
    utr_rect <- make_rect(tx_info[grep("utr", tx_info$feature), ], 
                          height=0.25, 
                          ann_ax)
    ncRNA_rect <- make_rect(tx_info[tx_info$feature == 'ncRNA', ],
                            height=0.25, 
                            ann_ax)
    intron_arrow <- make_arrows(tx_info[tx_info$feature == 'intron', ], ann_ax, 
                                arrowlen = width(target_range) * 0.01)
    # Compile new shapes
    tx_shapes <- c(cds_rect, utr_rect, ncRNA_rect, intron_arrow)
    
    # Add shapes to list if it already exists
    if(!is.null(plotly_obj$x$layout$shapes)){
      tx_shapes <- c(plotly_obj$x$layout$shapes, tx_shapes)
    }
    plotly_obj$x$layout$shapes <- tx_shapes
  }
  return(plotly_obj)
}

#' make_rect
#' Make the rectangle objects. These are used to display annotation features such
#' as cds, and utr.
#'
#' @param df data.frame: data used to draw the shapes. Required columns include 
#' "start", "end", and "stepping"
#' @param height numeric: the height of the rectangles
#' @param yref character: axis on which to draw the rectangles, in the form of 
#' "yaxis[2-9]" or "y[2-9]"
#' @param fillcolor character: color of the rectangle
#'
#' @return
#'
#' @author Justin Finkle
#' @examples
make_rect <- function(df, 
                      height, 
                      yref,
                      fillcolor = 'lightslateblue'){
  fillcolor <- make_plotly_color(fillcolor)
  yref <- gsub("yaxis", "y", yref)
  if(nrow(df)>0){
    rect_list <- vector("list", nrow(df))
    for(e in 1:nrow(df)){
      row <- df[e, ]
      rect_list[[e]] <- list(type = "rect", fillcolor = fillcolor, opacity = 1, 
                             line=list(width=0),
                             x0 = row$start, x1 = row$end, xref = "x",
                             y0 = row$stepping-height, y1 = row$stepping+height, 
                             yref = yref)
    }
  } else {
    rect_list <- NULL
  }
  return(rect_list)
}

#' make_arrows
#' Draw arrows for introns
#'
#' @param df data.frame: data used to draw the shapes. Required columns include 
#' "start", "end", "strand", "midpoint", and "stepping"
#' @param yref character: axis on which to draw the rectangles, in the form of 
#' "yaxis[2-9]" or "y[2-9]"
#' @param arrowlen numeric: how long the arrow should be, in x coordinates
#' @param arrowheight numeric: how tall the arrow should be, in yref coordinates
#' @param arrowgap  numeric: gap between arrows on long introns
#'
#' @return
#'
#' @examples
make_arrows <- function(df, 
                        yref, 
                        arrowlen = 500, 
                        arrowheight = 0.15, 
                        arrowgap = 1500){
  # Set the y ref
  yref <- gsub("yaxis", "y", yref)
  if(nrow(df)>0){
    line_list <- vector("list", nrow(df))
    arrow_list <- vector("list", nrow(df))
    for(i in 1:nrow(df)){
      row <- df[i, ]
      
      # Make the intron line
      line_list[[i]] <- list(x0 = row$start, y0=row$stepping, 
                             x1 = row$end, y1 = row$stepping, 
                             xref = "x", yref = yref,
                             type = "line",
                             line = list(width = 0.5, 
                                         dash = ifelse(row$strand == "-", 
                                                       "dot","solid")))
      # Add arrows to the lines
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

#' arrow_helper
#' Draw lines that make arrows on the intron shapes
#'
#' @param arrow_start 
#' @param strand 
#' @param arrowlen 
#' @param arrowheight 
#' @param y 
#' @param yref 
#'
#' @return
#'
#' @examples
arrow_helper <- function(arrow_start, 
                         strand, arrowlen, 
                         arrowheight, 
                         y, 
                         yref){
  arrow_end <- ifelse(strand =="-", arrow_start + arrowlen, 
                      arrow_start - arrowlen)
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



# =============================================================================
# =============================================================================
# Deprecated Functions
# =============================================================================
# =============================================================================

# Draw arrows as annotations. This has unexpected behavior, and drawing as 
# shapes is preferred.
old_make_arrows <- function(df, yref){
  yref <- gsub("yaxis", "y", yref)
  if(nrow(df)>0){
    arrow_list <- vector("list", nrow(df))
    for(i in 1:nrow(df)){
      row <- df[i, ]
      xstart <- ifelse(row$strand =="-", row$start, row$end)
      xend <- ifelse(row$strand =="-", row$end, row$start)
      arrow_list[[i]] <- list(x = xstart, y=row$stepping, xref = "x", 
                              yref = yref,
                              showarrow = TRUE, ax = xend, ay=row$stepping,
                              axref='x', ayref= yref, arrowwidth = 1, text="")
      
    }
  } else {
    arrow_list <- NULL
  }
  return(arrow_list)
}


# Crop introns so shape only falls within view range. Not needed now that intron
# arrows are drawn as shapes, not annotations
crop_introns <- function(introns, target_range){
  introns$start <- pmax(introns$start, start(target_range))
  introns$end <- pmin(introns$end, end(target_range))
  introns$midpoint <- pmin(pmax(introns$midpoint, 
                                start(target_range)), end(target_range))
  return(introns)
}

coverage_heatmap <- function(cvg){
  hm_vals <- apply(t(as.matrix((mcols(cvg)))), 2, rev)
  midpoint <- get_midpoint(cvg)
  p <- plot_ly(z=hm_vals, y=rev(colnames(mcols(cvg))), x=midpoint, 
               type='heatmap',
               hoverinfo='x+z', 
               colorscale = continuous_colorscale("Purples")(hm_vals))
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
  mcols(dense_rna)$feature <- 
    splitAsList(mcols(rna)$feature[queryHits(overlaps)],
                factor(subjectHits(overlaps)))
  mcols(dense_rna)$transcript <- 
    splitAsList(mcols(rna)$transcript[queryHits(overlaps)],
                factor(subjectHits(overlaps)))
  return(c(dense_introns, dense_rna))
}