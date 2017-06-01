info_dir <- "/gne/research/workspace/rinaldj1/projects/bulkChipPipeline/data/CARROLL_PR/data/"
s_file <- paste(info_dir, "samples.txt", sep='')
p_file <- paste(info_dir, "pairs.txt", sep='')
fi <- chipVis::get_chip_file_info(s_file, p_file)
org <- "human"
genome <- "GRCh38"

anno <- readr::read_csv("/gne/research/workspace/rinaldj1/projects/bulkChipPipeline/data/CARROLL_PR/datDex.csv")[,c("run","sample_attribute")]

sample_attributes = lapply(strsplit(anno$sample_attribute, split = "||", fixed = TRUE),
                           function(x) setNames(lapply(x, function(y) stringr::str_trim(strsplit(y,":")[[1]][2])),
                                                stringr::str_trim(sapply(x, function(y) strsplit(y,":")[[1]][1]))))

sample_attributes_df <- plyr::ldply(sample_attributes, data.frame, stringsAsFactors = FALSE)

fi_attributes <- as.data.frame(t(sapply(fi$Sample.Name, function(x){
  if (x %in% anno$run){
    return(sample_attributes_df[which(anno$run == x), ])
  } else{
    return(sample_attributes_df[grep(x,sample_attributes_df$source_name)[1], ])
  }
})))

fi <- cbind(fi, fi_attributes)

saveRDS(fi,
        file = "/gne/research/workspace/schepa/Carroll_Data/chip_file_info.Rds")



rpkmEset <- readRDS("/gne/research/workspace/schepa/Carroll_Data/GSE68358.rpkmEset.Rds")

tr <- transcriptsBy(TxDb.Hsapiens.BioMart.igis, by = "gene") %>% unlist()
tr$entrez <- stringr::str_extract(names(tr),"(?<=:)(.*)")
tr$symbol <- annotate::getSYMBOL(tr$entrez, data='org.Hs.eg.db')
names(tr) = NULL
trdf <- as.data.frame(tr)

truniq <- trdf %>% group_by(entrez) %>% summarize(chr = first(seqnames),
                                                  start = min(start),
                                                  end = min(end),
                                                  strand = first(strand),
                                                  symbol = first(symbol))

fData(rpkmEset) <- truniq

saveRDS(rpkmEset,
        file = "/gne/research/workspace/schepa/Carroll_Data/GSE68358.rpkmEset.with_fdata.Rds")



sf <- get_scaling_factor(fi)[2,]







