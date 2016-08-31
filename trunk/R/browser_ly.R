coverage_heatmap <- function(cvg){
  hm_vals <- apply(t(as.matrix((mcols(cvg)))), 2, rev)
  midpoint <- get_midpoint(cvg)
  p <- plot_ly(z=hm_vals, y=rev(colnames(mcols(cvg))), x=midpoint, type='heatmap',
               hoverinfo='x+z', colorscale = continuous_colorscale("Purples")(hm_vals))
  return(p)
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
    b_plot <- plot_browserly_tracks(target_range = target_range, 
                                    tx_data = tx_data,
                                    cvg = cvg_gr,
                                    hm_thresh = hm_thresh,
                                    type = type)
    return(b_plot)
  }
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
plot_browserly_tracks <- function(target_range, tx_data, cvg,
                        hm_thresh = 4,
                        stacking = 'dense',
                        type = NULL){
  # Make GeneRegion plot
  tx_info <- get_tx_annotation(range=target_range, tx_data = tx_data)
  
  # Prepare tx info for plotting
  tx_info <- set_tx_level(tx_info, stacking = stacking)
  tx_info <- biovizBase::mold(tx_info)
  
  # Make "invisible" transcript plot
  tx_track <- browserly_tx_track(tx_info = tx_info)
  
  plots <- list("GeneRegion"= tx_track)
  # Make supblots
  # Decide the type of datatrack to plot if not provided
  if(is.null(type)){
    type <- ifelse(ncol(mcols(cvg)) > hm_thresh, "heatmap", "hist")
  }
  
  cvg_track <- browserly_cvg_track(cvg = cvg, target_range = target_range,
                                   type = type)
  
  plots <- modifyList(plots, cvg_track)
  track_heights <- c(0.3, rep(0.7/length(cvg_track), length(cvg_track)))
  sp <- subplot(plots, nrows=length(plots), shareX = TRUE, heights=track_heights)
  
  # Modify layouts
  trace_names <- sapply(sp$x$data, function(x){x$name})
  trace_axes <- lapply(sp$x$data, function(x) {gsub("y", "yaxis", x$yaxis)})
  names(trace_axes) <- trace_names                                                                                                     
  grt_ax <- trace_axes[["GeneRegion"]]
  
  letterwidth <- 9
  maxlabellength <- max(sapply(colnames(mcols(cvg)), nchar))
  if(type == 'heatmap'){
    left_margin = max(letterwidth*maxlabellength, 60)
    margin = list(l=left_margin)
  } else {
    margin =  list()
  }
  sp <- add_tx_shapes(sp, tx_info, target_range, grt_ax)
  sp <- modify_y(sp, trace_axes, grt_ax, type)
  sp <- sp %>% layout(xaxis=list(range=c(start(target_range), end(target_range))),
                      margin=margin)
  return(sp)
  
}

modify_y <- function(plotly_obj, trace_axes, grt_ax, type){
  # Modify transcript axis
  plotly_obj$x$layout[[grt_ax]][['ticks']] <- ""
  plotly_obj$x$layout[[grt_ax]][['showticklabels']] <- FALSE
  plotly_obj$x$layout[[grt_ax]][['title']] <- "Transcripts"
  plotly_obj$x$layout[[grt_ax]][['showgrid']] <- FALSE
  
  data_axes <- trace_axes[trace_axes!=grt_ax]
  
  # Scale the data for the histograms
  if(type!='heatmap'){
    score_max <- max(unlist(lapply(plotly_obj$x$data, function(x) 
      {if(class(x$y)=='numeric') return(x$y)})))
  }
  for(sample in names(data_axes)){
    ax <- data_axes[[sample]]
    if(type!='heatmap'){
      plotly_obj$x$layout[[ax]][['range']] <- c(0, score_max)
      plotly_obj$x$layout[[ax]][['title']] <- sample
      
      # plotly_obj$x$layout[[ax]][['titlefont']]2 <- list(size=10)
    } else {
      plotly_obj$x$layout[[ax]]['ticks'] <- ""
    }
  }
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

browserly_tx_track <- function(tx_info, 
                               hover_info = c('cds', 'utr3', 'utr5', 'ncRNA')){
  
  # Initizialize Plot
  tx_track <- plot_ly(type='scatter', name = 'GeneRegion')
  for(h_feature in hover_info){
    h <- tx_info[tx_info$feature == h_feature, ]
    tx_track <- tx_track %>% add_trace(x=h$midpoint, y=h$stepping, opacity=0,
                                       text = h$transcript, 
                                       name=h_feature,
                                       hoverinfo="x+name+text", showlegend=FALSE)
  }
  
  return(tx_track)
}

browserly_cvg_track <- function(cvg, target_range, 
                                type = 'heatmap',
                                colors = NULL){
  
  cvg_track <- list()
  if(type == 'heatmap'){
    cvg_track[['heatmap']] <- coverage_heatmap(cvg)
  }
  else{
    midpoint <- get_midpoint(cvg)
    for(n in colnames(mcols(cvg))){
      cvg_track[[n]] <- plot_ly(x=midpoint, y=mcols(cvg)[,n], type='scatter',
                                fill='tozeroy', mode='', showlegend=FALSE, name=n,
                                hoverinfo='x+y')
    }
  }
  return(cvg_track)
}

set_tx_level <- function(tx_gr, 
                         stacking = c('dense', 'squish')){
  stacking <- match.arg(stacking)
  if(length(tx_gr)){
    if(stacking == 'squish'){
      tx_gr <- biovizBase::addStepping(tx_gr, group.name = "transcript",
                                       group.selfish = FALSE)
    } else if(stacking == 'dense'){
      # tx_gr <- collapse_tx(tx_gr)
      mcols(tx_gr)$stepping <- 1
    }
    return(tx_gr)
  }
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
  
  
