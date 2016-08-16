coverage_heatmap <- function(b_cvg){
  hm_vals <- apply(t(as.matrix(as.data.frame(mcols(b_cvg)))), 2, rev)
  bin_df <- biovizBase::mold(b_cvg)
  p <- plot_ly(z=hm_vals, y=rownames(hm_vals), x=bin_df$midpoint, type='heatmap',
               hoverinfo='x+z', colorscale = continuous_colorscale("Purples")(hm_vals))
  return(p)
}

crop_introns <- function(introns, range){
  introns$start <- pmax(introns$start, start(range))
  introns$end <- pmin(introns$end, end(range))
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


#' Title
#'
#' @param cvg_files 
#' @param genome 
#' @param tx_data 
#' @param hm_thresh 
#' @param type 
#' @param binsize 
#' @param cvg_scaling 
#'
#' @return
#' @export
#'
#' @examples
make_browserly_function <- function(cvg_files, 
                                    tx_data,
                                    hm_thresh = 4,
                                    type = NULL,
                                    binsize = 1000,
                                    cvg_scaling = NULL){
  plot_browserly <- function(target_range){
    cvg_list <- get_coverage_in_range(bwList = cvg_files,
                                      target_range = target_range, 
                                      names = names(cvg_files), 
                                      cvg_scaling = cvg_scaling)
    b_plot <- plot_browserly_tracks(range = target_range, 
                                    tx_data = tx_data,
                                    cvg_list = cvg_list,
                                    hm_thresh = hm_thresh,
                                    type = type,
                                    binsize = binsize)
    return(b_plot)
  }
  return(plot_browserly)
}


#' Title
#'
#' @param txdb 
#' @param range 
#' @param tx_data 
#' @param cvg_list 
#' @param hm_thresh 
#' @param type 
#' @param binsize 
#'
#' @return
#' @export
#'
#' @examples
plot_browserly_tracks <- function(range, tx_data, cvg_list,
                        hm_thresh = 4,
                        type = NULL,
                        binsize = 1000){
  # Make GeneRegion plot
  tx_info <- get_tx_annotation(range=range, tx_data = tx_data)
  
  # Prepare tx info for plotting
  if(length(tx_info)){
    tx_info <- biovizBase::addStepping(tx_info, group.name = "transcript",
                                       group.selfish = FALSE)
  }
  tx_info <- biovizBase::mold(tx_info)
  
  # Make "invisible" transcript plot
  tx_track <- browserly_tx_track(tx_info = tx_info)
  
  plots <- list("GeneRegion"= tx_track)
  # Make supblots
  # Decide the type of datatrack to plot if not provided
  if(is.null(type)){
    type <- ifelse(length(cvg_list) > hm_thresh, "heatmap", "hist")
  }
  
  cvg_track <- browserly_cvg_track(cvg_list=cvg_list, range = range,
                                   type = type, binsize = binsize)
  
  plots <- modifyList(plots, cvg_track)
  track_heights <- c(0.3, rep(0.7/length(cvg_track), length(cvg_track)))
  sp <- subplot(plots, nrows=length(plots), shareX = TRUE, heights=track_heights)
  
  # Modify layouts
  trace_names <- sapply(sp$x$data, function(x){x$name})
  trace_axes <- lapply(sp$x$data, function(x) {gsub("y", "yaxis", x$yaxis)})
  names(trace_axes) <- trace_names                                                                                                     
  grt_ax <- trace_axes[["GeneRegion"]]
  
  letterwidth <- 9
  maxlabellength <- max(sapply(names(cvg_list), nchar))
  if(type == 'heatmap'){
    left_margin = max(letterwidth*maxlabellength, 60)
    margin = list(l=left_margin)
  } else {
    margin =  list()
  }
  sp <- add_tx_shapes(sp, tx_info, range, grt_ax)
  sp <- modify_y(sp, cvg_list, trace_axes, grt_ax, type)
  sp <- sp %>% layout(xaxis=list(range=c(start(range), end(range))),
                      margin=margin)
  return(sp)
  
}

modify_y <- function(plotly_obj, cvg_list, trace_axes, grt_ax, type){
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

add_tx_shapes <- function(plotly_obj, tx_info, range, grt_ax){
  # Add annotations
  cds_rect <- make_rect(tx_info[tx_info$feature == 'cds', ], height = 0.4, grt_ax)
  utr_rect <- make_rect(tx_info[grep("utr", tx_info$feature), ], height=0.2, grt_ax)
  ncRNA_rect <- make_rect(tx_info[tx_info$feature == 'ncRNA', ], height=0.2, grt_ax)
  
  # Crop arrows to view range
  cropped_introns <- crop_introns(tx_info[tx_info$feature == 'intron', ], range)
  intron_arrow <- make_arrows(cropped_introns, grt_ax)
  grt_layout <- list(showlegend=FALSE, shapes=c(cds_rect, utr_rect, ncRNA_rect),
                     annotations=intron_arrow)
  plotly_obj$x$layout <- modifyList(plotly_obj$x$layout, grt_layout)
  grt_y_layout <- list(autorange='reversed', showticklabels=FALSE, showticks=FALSE,
                       title='Transcripts')
  plotly_obj$x$layout[[grt_ax]] <- modifyList(plotly_obj$x$layout[[grt_ax]], grt_y_layout)
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

browserly_cvg_track <- function(cvg_list, range, 
                                binsize = 1000,
                                type = 'heatmap',
                                colors = NULL){
  bin_cvg <- bin_coverage_in_range(cvg_list, range, binwidth = binsize)
  bin_df <- biovizBase::mold(bin_cvg)
  cvg_track <- list()
  if(type == 'heatmap'){
    cvg_track[['heatmap']] <- coverage_heatmap(bin_cvg)
  }
  else{
    for(n in names(cvg_list)){
      cvg_track[[n]] <- plot_ly(x=bin_df$midpoint, y=bin_df[[n]], type='scatter',
                                fill='tozeroy', mode='', showlegend=FALSE, name=n,
                                hoverinfo='x+y')
    }
  }
  return(cvg_track)
}
