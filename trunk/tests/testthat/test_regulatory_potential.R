context('Regulatory Potential')

test_that("regulatory potential score works for GRanges tss input", {
  p <- GenomicRanges::GRanges("chr1", 
                             IRanges::IRanges(start = c(1,100, 500000), width = 500))
  t <- GenomicRanges::GRanges("chr1", 
                              IRanges::IRanges(start = c(5,501000, 800000), width = 1), 
                              strand = c("+","-","+"))
  rp <- compute_regulatory_potential(p, t, max_dist = 100000)
  tf <- function(x){
    sum(exp(-(0.5 + 4 * x / 100000)))
  }
  expect_equal(length(rp), length(t))
  expect_equal(rp[1], tf(c(0,94)))
  expect_equal(rp[2], tf(500))
  expect_equal(rp[3], 0)
})


test_that("regulatory potential score works for GRangesList tss input", {
  p <- GenomicRanges::GRanges("chr1", 
                              IRanges::IRanges(start = c(1,100, 500000), width = 500))
  t <- GenomicRanges::GRangesList(GenomicRanges::GRanges("chr1", 
                              IRanges::IRanges(start = c(5,10), width = 1), 
                              strand = c("+","+")),
                              GenomicRanges::GRanges("chr1", 
                                                     IRanges::IRanges(start = c(501000), width = 1), 
                                                     strand = c("-")))
  rp <- compute_regulatory_potential(p, t, max_dist = 100000)
  tf <- function(x){
    sum(exp(-(0.5 + 4 * x / 100000)))
  }
  expect_equal(length(rp), length(t))
  expect_equal(rp[1], tf(c(0,89)))
  expect_equal(rp[2], tf(500))
})

