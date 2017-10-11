context("iheatmap")

data("rpkm_chr21", package = "GenomicWidgets")
rpkm_sub <- rpkm_chr21[1:20,1:10]


test_that("iheatmap works with SummarizedExperiment and minimal args",{
  
  se_ihm <- iheatmap(rpkm_sub)
  
  expect_iheatmap(se_ihm, "expression_heatmap_basic")
  
})

test_that("iheatmap works with SummarizedExperiment with anno",{
  
  se_ihm <- iheatmap(rpkm_sub, "rpkm",
           x = colData(rpkm_sub)$STD_NAME, 
           y = rowData(rpkm_sub)$SYMBOL, 
           col_annotation = colData(rpkm_sub)[,c("TYPE","SEX")])
  
  expect_iheatmap(se_ihm, "expression_heatmap")
  
})

test_that("add_iheatmap works with SummarizedExperiment",{
  
  se_ihm <- iheatmap(rpkm_sub, "rpkm",
                     x = colData(rpkm_sub)$STD_NAME, 
                     y = rowData(rpkm_sub)$SYMBOL, 
                     col_annotation = colData(rpkm_sub)[,c("TYPE","SEX")])
  
  add_ihm <- add_iheatmap(se_ihm, rpkm_sub, "rpkm",
                          x = colData(rpkm_sub)$STD_NAME, 
                          col_annotation = colData(rpkm_sub)[,c("TYPE","SEX")])
  
  expect_iheatmap(add_ihm, "expression_heatmap_double")
  
})

test_that("add_iheatmap vertical works with SummarizedExperiment",{
  
  se_ihm <- iheatmap(rpkm_sub, "rpkm",
                     x = colData(rpkm_sub)$STD_NAME, 
                     y = rowData(rpkm_sub)$SYMBOL, 
                     col_annotation = colData(rpkm_sub)[,c("TYPE","SEX")],
                     orientation = "vertical")
  
  add_ihm <- add_iheatmap(se_ihm, rpkm_sub, "rpkm",
                          y = rowData(rpkm_sub)$SYMBOL)
  
  expect_iheatmap(add_ihm, "expression_heatmap_double_vertical","vertical")
  
})