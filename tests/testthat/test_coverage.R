context("Coverage")

test_that("make_coverage_tracks works on uniform binsize=100 for both methods a and j", {
  fi <- c('../testdata/test_coverage1.bw', '../testdata/test_coverage2.bw')
  target_range <- get_view_range(chr='chr1', start=1, end=900)
  
  cvg <- make_coverage_tracks(fi, target_range=target_range, binsize=100, method='a')
  expect_that(cvg$V1, equals(c(8,10,8271,1,6,6,386,39,20)))
  expect_that(cvg$V2, equals(c(5, 0.5, 4, 1, 1, 1, 1, 1, 3.6)))
  
  cvg <- make_coverage_tracks(fi, target_range=target_range, binsize=100, method='j')
  expect_that(unname(cvg$'1'), equals(c(8,10,8271,1,6,6,386,39,20)))
  expect_that(unname(cvg$'2'), equals(c(5, 0.5, 4, 1, 1, 1, 1, 1, 3.6)))
})

test_that("make_coverage_tracks works on unevenly distributed binsizes for both methods a and j", {
  fi <- c('../testdata/test_coverage1.bw', '../testdata/test_coverage2.bw')
  target_range <- get_view_range(chr='chr1', start=1, end=850)
  
  cvg_a <- make_coverage_tracks(fi, target_range=target_range, binsize=100, method='a')
  expect_that(cvg_a$V1[1], equals(8))
  expect_that(cvg_a$V1[2], equals(9.87234, tol=0.0001))
  expect_that(cvg_a$V2[3], equals((2*33+6*33+(283-266)*4)/(283-188)))
  expect_that(cvg_a$V2[4], equals((4*(300-283)+1*(377-300))/(377-283)))
  
  cvg_j <- make_coverage_tracks(fi, target_range=target_range, binsize=100, method='j')
  expect_that(cvg_a$V1, is_equivalent_to(cvg_j$'1'))
  expect_that(cvg_a$V2, is_equivalent_to(cvg_j$'2'))
  })

test_that("make_coverage_tracks creates relatively uniform binsizes", {
  fi <- c('../testdata/test_coverage1.bw', '../testdata/test_coverage2.bw')
  target_range <- get_view_range(chr='chr1', start=1, end=860)
  
  cvg_j <- make_coverage_tracks(fi, target_range=target_range, binsize=100, method='j')
  expect_that(width(ranges(cvg_j)), equals(rep(94, 9), tolerance=1))
  
  cvg_a <- make_coverage_tracks(fi, target_range=target_range, binsize=100, method='a')
  expect_that(width(ranges(cvg_a)), equals(rep(94, 9), tolerance=1))
              
  expect_that(cvg_a$V1, is_equivalent_to(cvg_j$'1'))
  expect_that(cvg_a$V2, is_equivalent_to(cvg_j$'2'))
})