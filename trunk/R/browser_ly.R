txdb_to_shapes <- function(txdb, ranges){
  stopifnot(length(ranges) == 1)
  tx <- GenomicFeatures::transcriptsByOverlaps(txdb, ranges)
  expanded_range <- GenomicRanges::GRanges(seqnames(ranges), IRanges::IRanges(min(start(tx)), max(end(tx))))
  
  gr <- biovizBase::crunch(txdb, which = expanded_range)
  
  gr <- biovizBase::addStepping(gr, group.name = "tx_id", group.selfish = FALSE)
  
  df <- biovizBase::mold(gr)
  
  out <- list(exons = df[which(df$type == "exon"),], 
              introns = df[which(df$type == "gap"),],
              utr = df[which(df$type == "utr"),])
  return(out)
}

make_exons <- function(df){
  rect_list <- vector("list", nrow(df))
  for(e in 1:nrow(df)){
    row <- df[e, ]
    rect_list[[e]] <- list(type = "rect", fillcolor = "blue", opacity = 1, line=list(width=0),
                         x0 = row$start, x1 = row$end, xref = "x",
                         y0 = row$stepping-0.4, y1 = row$stepping+0.4, yref = "y")
  }
  return(rect_list)
}

make_introns <- function(df){
  arrow_list <- vector("list", nrow(df))
  for(i in 1:nrow(df)){
    row <- df[i, ]
    xstart <- ifelse(row$strand =="-", row$start, row$end)
    xend <- ifelse(row$strand =="-", row$end, row$start)
    arrow_list[[i]] <- list(x = xstart, y=row$stepping, xref = "x", yref = "y",
                            showarrow = TRUE, ax = xend, ay=row$stepping,
                            axref='x', ayref='y', arrowwidth =1)
    
  }
  return(arrow_list)
}

make_utr <- function(df){
  rect_list <- vector("list", nrow(df))
  for(e in 1:nrow(df)){
    row <- df[e, ]
    rect_list[[e]] <- list(type = "rect", fillcolor = "red", opacity = 1, line=list(width=0),
                           x0 = row$start, x1 = row$end, xref = "x",
                           y0 = row$stepping-0.2, y1 = row$stepping+0.2, yref = "y")
  }
  return(rect_list)
}

plot_gr <- function(tx_list, range){
  
  p <- plot_ly(tx_list$exons, x=~midpoint, y=~stepping, type='scatter', opacity=0,
               text = ~tx_name)
  exons <- make_exons(tx_list$exons)
  introns <- make_introns(tx_list$introns)
  utr <- make_utr(tx_list$utr)
  
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
    showline = TRUE,
    showticklabels = FALSE,
    showgrid = FALSE,
    autorange = 'reversed'
  )
  
  p <- p %>% layout(shapes = c(exons, utr), annotations=introns,
                    xaxis = ax, yaxis=ay)
  return(p)
}