context("tracks")

library(TxDb.Hsapiens.UCSC.hg19.knownGene)

track_plotter <- make_track_plotter(samp.info$fileName[1:3], 
                                    annotation = 
                                      TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                    track_names = samp.info$sampleName[1:3], 
                                    share_y = TRUE)

test_that("track plot generator is function",{
  
  expect_is(track_plotter,"function")

})

test_that("track plot generator works with single locus",{
  
  p1 <- track_plotter(resize(ctcf.peaks[1], width = 5000, fix = "center"))
  
  expect_is(p1, "LocusViewList")
  expect_genomic_widget(p1,"track_single_locus")
  
})

test_that("track plot generator works with multiple loci",{
  
  p1 <- track_plotter(resize(ctcf.peaks[1:3], width = 5000, fix = "center"))
  
  expect_is(p1, "LocusViewList")
  expect_genomic_widget(p1,"track_multi_locus")
  
})