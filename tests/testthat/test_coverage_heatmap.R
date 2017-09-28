context("coverage_heatmap")

test_mats <- readRDS("test_mats.Rds")

test_that("Coverage heatmap works with one assay from SummarizedExperiment",{
  
  cov_hm <- coverage_heatmap(test_mats,"Ctcf")
  
  expect_iheatmap(cov_hm, "coverage_single_se")
  
})

test_that("Coverage heatmap works with multiple assays from SummarizedExperiment",{
  
  cov_hm <- coverage_heatmap(test_mats,c("Ctcf","P300"))
  
  expect_iheatmap(cov_hm, "coverage_two_se")
  
})

test_that("Coverage heatmap works with matrix",{
  
  cov_hm <- coverage_heatmap(assays(test_mats)$Ctcf)
  
  expect_iheatmap(cov_hm, "coverage_single_matrix")
  
})

test_that("Coverage heatmap works with list of matrices",{
  
  cov_hm <- coverage_heatmap(as.list(assays(test_mats))[1:3])
  
  expect_iheatmap(cov_hm, "coverage_matrix_list")
  
})


geno_mats <- genomation::ScoreMatrixList(targets = samp.info$fileName[1:3],
                                         windows = resize(ctcf.peaks[1:10], 
                                                          width = 500, 
                                                          fix = "start"), 
                                         bin.num = 20)

test_that("Coverage heatmap works with ScoreMatrix",{
  
  cov_hm <- coverage_heatmap(geno_mats[[1]], -250, 250)
  
  expect_iheatmap(cov_hm, "coverage_ScoreMatrix")
  
})

test_that("Coverage heatmap works with ScoreMatrixList",{
  
  cov_hm <- coverage_heatmap(geno_mats, -250, 250)
  
  expect_iheatmap(cov_hm, "coverage_ScoreMatrixList")
  
})

## Add coverage heatmap

test_that("Adding Coverage heatmap works with one assay from SummarizedExperiment",{
  
  cov_hm <- coverage_heatmap(test_mats,"Ctcf")
  
  add_hm <- add_coverage_heatmap(cov_hm, test_mats, "P300")
  
  expect_iheatmap(add_hm, "add_coverage_single_se")
  
})

test_that("Adding Coverage heatmap works with 2 assays from SummarizedExperiment",{
  
  cov_hm <- coverage_heatmap(test_mats,"Ctcf")
  
  add_hm <- add_coverage_heatmap(cov_hm, test_mats, c("P300","Znf143"))
  
  expect_iheatmap(add_hm, "add_coverage_two_se")
  
})


test_that("Adding Coverage heatmap works with matrix",{
  
  cov_hm <- coverage_heatmap(test_mats,"Ctcf")
  
  add_hm <- add_coverage_heatmap(cov_hm, assays(test_mats)$P300)
  
  expect_iheatmap(add_hm, "add_coverage_matrix")
  
})

test_that("Adding Coverage heatmap works with list of matrices",{
  
  cov_hm <- coverage_heatmap(test_mats,"Ctcf")
  
  add_hm <- add_coverage_heatmap(cov_hm, as.list(assays(test_mats)[c("P300","Znf143")]))
  
  expect_iheatmap(add_hm, "add_coverage_matrix_list")
  
})


test_that("Adding Coverage heatmap works with ScoreMatrix",{
  
  cov_hm <- coverage_heatmap(test_mats,"Ctcf")
  
  add_hm <- add_coverage_heatmap(cov_hm, geno_mats[[2]], -250, 250)
  
  expect_iheatmap(add_hm, "add_coverage_ScoreMatrix")
  
})

test_that("Adding Coverage heatmap works with ScoreMatrixList",{
  
  cov_hm <- coverage_heatmap(test_mats,"Ctcf")
  
  add_hm <- add_coverage_heatmap(cov_hm, geno_mats[2:3], -250, 250)
  
  expect_iheatmap(add_hm, "add_coverage_ScoreMatrixList")
  
})
