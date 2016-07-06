context('File Access')

test_that("loadNarrowPeaks returns GRanges", {
  gr_narrowpeak <- loadNarrowPeaks("../testdata/test.narrowPeak")
  expect_is(gr_narrowpeak, "GRanges")
})