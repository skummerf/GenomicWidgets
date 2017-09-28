# Some utility functions 

# makes x based on colnames of mat if available
# if not available, just uses 1 to number of columns
default_x <- function(mat){
  if (is.null(colnames(mat))){
    return(seq_len(ncol(mat)))
  } else{
    colnames(mat)
  }
}
# makes y based on rownames of mat if available
# if not available, just uses 1 to number of rows
default_y <- function(mat){
  if (is.null(rownames(mat))){
    return(seq_len(nrow(mat)))
  } else{
    rownames(mat)
  }
}
