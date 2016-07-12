context('File Access')

test_that("getChipFileInfo loads proper dataframe", {
  samples <- "../testdata/test_samples.txt"
  pairs <- "../testdata/test_pairs.txt"
  file_info <- getChipFileInfo(samples, pairs)
  expect_equal_to_reference(file_info, "../testdata/test_file_info.rds")
  rm(samples, pairs, file_info)
})