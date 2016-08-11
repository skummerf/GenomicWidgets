
trf <- make_track_function(fi$File.bw[unlist(T47D[1:3])],
                           genome, 
                           TxDb.Hsapiens.BioMart.igis, 
                           txd, 
                           sf[unlist(T47D[1:3])])

trf2 <- make_browserly_function(fi$File.bw[unlist(T47D[1:3])],
                            TxDb.Hsapiens.BioMart.igis, 
                            txd, 
                            cvg_scaling = sf[unlist(T47D[1:3])])

hmc <- heatmap_click(rna_rp_hm, tss_ranges[sig])


heatmap_to_browserly_shiny(rna_rp_hm, trf2, hmc)

