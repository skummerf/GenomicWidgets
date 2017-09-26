context("Coverage matrix")

test_that("Coverage matrix works with bam file",{
  test_mats <- make_coverage_matrix(samp.info$fileName[1:5], 
                                    ctcf.peaks[1:10], 
                                    input_names = samp.info$sampleName[1:5],
                                    up = 250, 
                                    down = 250, 
                                    binsize = 25)
  expect_is(test_mats, "SummarizedExperiment")
})