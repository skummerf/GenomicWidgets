context('File Access')

test_that("get_chip_file_info loads proper dataframe", {
  samples <- "../testdata/test_samples.txt"
  pairs <- "../testdata/test_pairs.txt"
  file_info <- get_chip_file_info(samples, pairs)
  expect_equal_to_reference(file_info, "../testdata/test_file_info.rds")
  rm(samples, pairs, file_info)
})