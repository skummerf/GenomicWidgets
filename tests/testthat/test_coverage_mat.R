context("Coverage matrix")


geno_mats <- genomation::ScoreMatrixList(targets = samp.info$fileName[1:3],
                                           windows = resize(ctcf.peaks[1:10], 
                                                            width = 500, 
                                                            fix = "start"), 
                                           bin.num = 20)
test_that("Coverage matrix works with bam file",{
  
  test_mats <- make_coverage_matrix(samp.info$fileName[1:3], 
                                    ctcf.peaks[1:10], 
                                    input_names = samp.info$sampleName[1:3],
                                    up = 250, 
                                    down = 250, 
                                    binsize = 25)
  
  expect_is(test_mats, "SummarizedExperiment")
  
  expect_equivalent(assays(test_mats)[[1]], signif(geno_mats[[1]]@.Data,3))
  expect_equivalent(assays(test_mats)[[2]], signif(geno_mats[[2]]@.Data,3))
  expect_equivalent(assays(test_mats)[[3]], signif(geno_mats[[3]]@.Data,3))
  
})

test_that("Coverage matrix works with bigwig file",{
  
  test_mats <- make_coverage_matrix(c("bw1.bw","bw2.bw","bw3.bw"), 
                                    ctcf.peaks[1:10], 
                                    input_names = samp.info$sampleName[1:3],
                                    up = 250, 
                                    down = 250, 
                                    binsize = 25)
  
  expect_is(test_mats, "SummarizedExperiment")
  
  expect_equivalent(assays(test_mats)[[1]], signif(geno_mats[[1]]@.Data,3))
  expect_equivalent(assays(test_mats)[[2]], signif(geno_mats[[2]]@.Data,3))
  expect_equivalent(assays(test_mats)[[3]], signif(geno_mats[[3]]@.Data,3))
  
})


