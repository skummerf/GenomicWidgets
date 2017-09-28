context("iheatmap")

data(rpkm_chr21)

test_that("iheatmap works with SummarizedExperiment and minimal args",{
  
  se_ihm <- iheatmap(rpkm_chr21)
  
  expect_iheatmap(se_ihm, "expression_heatmap_basic")
  
})

test_that("iheatmap works with SummarizedExperiment with anno",{
  
  se_ihm <- iheatmap(rpkm_chr21, "rpkm",
           x = colData(rpkm_chr21)$STD_NAME, 
           y = rowData(rpkm_chr21)$SYMBOL, 
           col_annotation = colData(rpkm_chr21)[,c("TYPE","SEX")])
  
  expect_iheatmap(se_ihm, "expression_heatmap")
  
})

test_that("add_iheatmap works with SummarizedExperiment",{
  
  se_ihm <- iheatmap(rpkm_chr21, "rpkm",
                     x = colData(rpkm_chr21)$STD_NAME, 
                     y = rowData(rpkm_chr21)$SYMBOL, 
                     col_annotation = colData(rpkm_chr21)[,c("TYPE","SEX")])
  
  add_ihm <- add_iheatmap(se_ihm, rpkm_chr21, "rpkm",
                          x = colData(rpkm_chr21)$STD_NAME, 
                          col_annotation = colData(rpkm_chr21)[,c("TYPE","SEX")])
  
  expect_iheatmap(add_ihm, "expression_heatmap_double")
  
})

test_that("add_iheatmap vertical works with SummarizedExperiment",{
  
  se_ihm <- iheatmap(rpkm_chr21, "rpkm",
                     x = colData(rpkm_chr21)$STD_NAME, 
                     y = rowData(rpkm_chr21)$SYMBOL, 
                     col_annotation = colData(rpkm_chr21)[,c("TYPE","SEX")],
                     orientation = "vertical")
  
  add_ihm <- add_iheatmap(se_ihm, rpkm_chr21, "rpkm",
                          y = rowData(rpkm_chr21)$SYMBOL)
  
  expect_iheatmap(add_ihm, "expression_heatmap_double_vertical","vertical")
  
})