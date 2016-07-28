plot_tx <- function(txdb, range, tx_data){
  tx <- transcriptsByOverlaps(txdb, range)
  tx_names <- tx$tx_name
  gr <- get_plot_ranges(tx_names, tx_data)
  if(!length(gr)){ return(NULL)}
  gr <- biovizBase::addStepping(gr, group.name = "transcript", group.selfish = FALSE)
  gr <- biovizBase::mold(gr)
  
  cds <- gr[gr$feature == 'cds', ]
  p <- plot_ly(x=cds$midpoint, y=cds$stepping, type='scatter', opacity=0,
               text = cds$transcript)
  cds_rect <- make_rect(cds, height = 0.4)
  utr_rect <- make_rect(gr[grep("utr", gr$feature), ], height=0.2)
  intron_arrow <- make_arrows(gr[gr$feature == 'intron', ])
  
  # Display defaults
  ax <- list(
    title = "",
    showgrid = FALSE,
    range= c(start(range), end(range)),
    showline = TRUE,
    side = 'top'
  )
  
  ay <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE,
    autorange = 'reversed'
  )
  
  p <- p %>% layout(shapes = c(cds_rect, utr_rect), annotations=intron_arrow,
                    xaxis = ax, yaxis=ay)
  return(p)
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

make_rect <- function(df, height){
  rect_list <- vector("list", nrow(df))
  for(e in 1:nrow(df)){
    row <- df[e, ]
    rect_list[[e]] <- list(type = "rect", fillcolor = "blue", opacity = 1, line=list(width=0),
                         x0 = row$start, x1 = row$end, xref = "x",
                         y0 = row$stepping-height, y1 = row$stepping+height, yref = "y")
  }
  return(rect_list)
}

make_arrows <- function(df){
  if(nrow(df)>0){
    arrow_list <- vector("list", nrow(df))
    for(i in 1:nrow(df)){
      row <- df[i, ]
      xstart <- ifelse(row$strand =="-", row$start, row$end)
      xend <- ifelse(row$strand =="-", row$end, row$start)
      arrow_list[[i]] <- list(x = xstart, y=row$stepping, xref = "x", yref = "y",
                              showarrow = TRUE, ax = xend, ay=row$stepping,
                              axref='x', ayref='y', arrowwidth =1)
    
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
