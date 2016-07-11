context('File Access')

test_that("getChipFileInfo loads proper dataframe", {
  samples <- "../testdata/test_samples.txt"
  pairs <- "../testdata/test_paris.txt"
  file_info <- getChipFileInfo(samples, pairs)
  expect_is(file_info, "data.frame"
  rm(samples, pairs)
})