get_tx_annotation <- function(txdb, range, tx_data){
  tx <- transcriptsByOverlaps(txdb, range)
  tx_names <- tx$tx_name
  gr <- get_plot_ranges(tx_names, tx_data)
  if(length(gr)){
    gr <- biovizBase::addStepping(gr, group.name = "transcript", 
                                  group.selfish = FALSE)
  }
  df <- biovizBase::mold(gr)
  return(df)
}

get_plot_ranges <- function(tx_names, tx_data){
  if(!length(tx_names)){ return(GRanges())}
  for(n in names(tx_data)){
    if(length(tx_data[[n]])){
      tx_subset <- tx_names[tx_names %in% names(tx_data[[n]])]
      part_gr <- unlist(tx_data[[n]][tx_subset])
      part_gr$transcript <- names(part_gr)
      part_gr$feature <- n
      if(!('exon_name' %in% colnames(mcols(part_gr)))){
        part_gr$exon_name <- NA
      }
      part_gr <- part_gr[, c('transcript', 'feature', 'exon_name')]
    } else {
      part_gr <- GRanges()
    }
    if(exists("parts_list")){
      parts_list <- c(parts_list, part_gr)
    } else {
      parts_list <- part_gr
    }
  }
  return(parts_list)
}

make_rect <- function(df, height, y_idx){
  if(nrow(df)>0){
    rect_list <- vector("list", nrow(df))
    yref = paste0("y", y_idx)
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

make_arrows <- function(df, y_idx){
  if(nrow(df)>0){
    yref = paste0("y", y_idx)
    arrow_list <- vector("list", nrow(df))
    for(i in 1:nrow(df)){
      row <- df[i, ]
      xstart <- ifelse(row$strand =="-", row$start, row$end)
      xend <- ifelse(row$strand =="-", row$end, row$start)
      arrow_list[[i]] <- list(x = xstart, y=row$stepping, xref = "x", yref = yref,
                              showarrow = TRUE, ax = xend, ay=row$stepping,
                              axref='x', ayref= yref, arrowwidth =1)
    
    }
  } else {
    arrow_list <- NULL
  }
  return(arrow_list)
}

load_tx_data <- function(txdb){
  message("Loading txdb data...")
  return(list(intron = intronsByTranscript(txdb, use.names=TRUE),
              utr5 = fiveUTRsByTranscript(txdb, use.names=TRUE),
              utr3 = threeUTRsByTranscript(txdb, use.names=TRUE),
              cds = cdsBy(txdb, by='tx', use.names=TRUE),
              exon = exonsBy(txdb, by='tx', use.names=TRUE))
  )
}

make_tracks <- function(txdb, range, tx_data, cvg_gr){
  # Make GeneRegion plot
  tx_info <- get_tx_annotation(txdb = txdb, range=range, tx_data = tx_data)
  # Make supblots
  bcvg <- bin_coverage_in_range(cvg_gr, range.gene, binwidth = 100)
  bin_df <- biovizBase::mold(bcvg)
  samples <- names(mcols(bcvg))
  plots <- list()
  
  cds <- tx_info[tx_info$feature == 'cds', ]
  plots[["GeneRegion"]] <- plot_ly(x=cds$midpoint, y=cds$stepping, type='scatter', opacity=0,
                                   text = cds$transcript, name="GeneRegion")
  for(n in samples){
    plots[[n]] <- plot_ly(x=bin_df$midpoint, y=bin_df[[n]],
                          type='bar', showlegend=FALSE, name=n)
  }
  
  track_heights <- c(0.2, rep(0.8/length(samples), length(samples)))
  sp <- subplot(plots, nrows=length(samples)+1, shareX = TRUE, heights=track_heights)
  track_names <- sapply(sp$x$data, function(x){x$name})
  grt_idx <- which(track_names == "GeneRegion")
  grt_idx <- ifelse(grt_idx == 1, '', grt_idx)
  
  # Make shapes
  cds_rect <- make_rect(cds, height = 0.4, grt_idx)
  utr_rect <- make_rect(tx_info[grep("utr", tx_info$feature), ], height=0.2, grt_idx)
  intron_arrow <- make_arrows(tx_info[tx_info$feature == 'intron', ], grt_idx)
  
  # Modify layouts
  grt_layout <- list(showlegend=FALSE, shapes=c(cds_rect, utr_rect),
                     annotations=intron_arrow)
  sp$x$layout <- modifyList(sp$x$layout, grt_layout)
  grt_y_layout <- list(autorange='reversed', showticklabels=FALSE, showticks=FALSE,
                       title='Transcripts')
  sp$x$layout[[paste0("yaxis",grt_idx)]] <- modifyList(sp$x$layout[[paste0("yaxis",grt_idx)]], grt_y_layout)
  sp <- sp %>% layout(xaxis=list(range=c(start(range), end(range))))
  return(sp)
    # Coverage and Heatmap
      # Titles
      # Colors
      # Y axes
      # Groups
    # Annotation
  
    # SNPs
    # Collapsing and TSS
  
}


