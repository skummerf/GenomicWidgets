context("Coverage matrix")

test_that("Coverage matrix works with bam file",{
  
  test_mats <- make_coverage_matrix(samp.info$fileName[1:3], 
                                    ctcf.peaks[1:10], 
                                    input_names = samp.info$sampleName[1:3],
                                    up = 250, 
                                    down = 250, 
                                    binsize = 25)
  
  geno_mats <- genomation::ScoreMatrixList(targets = samp.info$fileName[1:3],
                                           windows = resize(ctcf.peaks[1:10], 
                                                            width = 500, 
                                                            fix = "start"), 
                                           bin.num = 20)
  
  expect_is(test_mats, "SummarizedExperiment")
  
  expect_equivalent(assays(test_mats)[[1]], geno_mats[[1]]@.Data)
  expect_equivalent(assays(test_mats)[[2]], geno_mats[[2]]@.Data)
  expect_equivalent(assays(test_mats)[[3]], geno_mats[[3]]@.Data)
  
})

