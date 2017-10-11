context("tracks")

library(TxDb.Hsapiens.UCSC.hg19.knownGene)

track_params <- set_track_parameters(samp.info$fileName[1:3], 
                                    annotation = 
                                      TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                    track_names = samp.info$sampleName[1:3], 
                                    share_y = TRUE)

test_that("track plot params is right class",{
  
  expect_is(track_params,"TrackParameters")

})

test_that("track plot generator works with single locus",{
  
  p1 <- plot_tracks(resize(ctcf.peaks[1], width = 5000, fix = "center"),
                    track_params)
  
  expect_is(p1, "LocusViewList")
  expect_genomic_widget(p1,"track_single_locus")
  
})

test_that("track plot generator works with multiple loci",{
  
  p1 <- plot_tracks(resize(ctcf.peaks[1:3], width = 5000, fix = "center"),
                    track_params)
  
  expect_is(p1, "LocusViewList")
  expect_genomic_widget(p1,"track_multi_locus")
  
})

test_that("track plot generator works with alternative offset",{
  
  p1 <- plot_tracks(resize(ctcf.peaks[1:3], width = 5000, fix = "center"),
                    track_params, offset = 0)
  
  expect_is(p1, "LocusViewList")
  expect_genomic_widget(p1,"track_multi_locus_offset")
  
})

test_that("track plot generator works with alternative xtitle",{
  
  p1 <- plot_tracks(resize(ctcf.peaks[1:3], width = 5000, fix = "center"),
                    track_params, xtitle = "CTCF peaks")
  
  expect_is(p1, "LocusViewList")
  expect_genomic_widget(p1,"track_multi_locus_xtitle")
  
})