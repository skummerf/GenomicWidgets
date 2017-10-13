context("track_summaries")
 
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

data("rpkm_chr21",package = "GenomicWidgets")
tss <- promoters(rowRanges(rpkm_chr21), up = 1, down = 1)


summary_params <- 
  set_summary_parameters(rpkm_chr21,
                         groups = "GROUP",
                         ranges = resize(tss, width = 5000, fix = "center"))

track_params <- set_track_parameters(samp.info$fileName[1:3], 
                                     annotation = 
                                       TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                     track_names = samp.info$sampleName[1:3], 
                                     share_y = TRUE,
                                     summary = summary_params)

test_that("track plus summary plot generator works w/ single locus use char",{
  
  p1 <- plot_tracks(rownames(rpkm_chr21)[1],
                    track_params)
  
  expect_is(p1, "GenomeTrackWidget")
  expect_genomic_widget(p1,"track_summary_single_locus")
  
})

test_that("track plus summary plot generator works w/ >1 loci using char",{
  
  p1 <- plot_tracks(rownames(rpkm_chr21)[1:3],
                    track_params)
  
  expect_is(p1, "GenomeTrackWidget")
  expect_genomic_widget(p1,"track_summary_multi_locus")
  
})

test_that("track plus summary plot generator works w/ single locus use range",{
  
  p1 <- plot_tracks(resize(tss[1], width = 5000, fix = "center"),
                    track_params)
  
  expect_is(p1, "GenomeTrackWidget")
  expect_genomic_widget(p1,"track_summary_single_locus_range_input")
  
})

test_that("track plus summary plot generator works w/ >1 loci using range",{
  
  p1 <- plot_tracks(resize(tss[1:3], width = 5000, fix = "center"),
                    track_params)
  
  expect_is(p1, "GenomeTrackWidget")
  expect_genomic_widget(p1,"track_summary_multi_locus_range_input")
  
})


