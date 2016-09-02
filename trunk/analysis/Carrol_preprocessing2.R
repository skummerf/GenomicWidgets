library(dplyr)
library(limma)
library(edgeR)
library(Biobase)
library(GenomicFeatures)
library(TxDb.Hsapiens.BioMart.igis3.0)
library(org.Hs.eg.db)

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

sf <- get_scaling_factor(fi)[2,]

sample_info <- cbind(fi[,c("Sample.Name","File.peaks","File.bam","File.bw","source_name","condition","cell.line")],
                     data.frame(library.size = sf))


#saveRDS(sample_info,file = "/gne/research/workspace/schepa/Carroll_Data/example_sample_info.Rds")



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

T47D_eset <- rpkmEset[,grep("T47D",pData(rpkmEset)$celltype)]
pData(T47D_eset)$treatment <- stringr::str_trim(sapply(stringr::str_split(pData(T47D_eset)$treatment,":"), 
                                        function(x) x[[2]])) %>% gsub("+","_",., fixed= TRUE) %>% 
  gsub("_3hr","",.)


treat_groups <- pData(T47D_eset)$treatment
dge <- DGEList(counts = exprs(T47D_eset)) %>% calcNormFactors()
design <- model.matrix(~ 0 + treat_groups)
colnames(design) = gsub("treat_groups", "", colnames(design))
v <- voom(dge, design=design)
fit <- lmFit(v, design)
contrasts <- makeContrasts(E2_Progesterone-E2, 
                           E2_R5020-E2, 
                           E2_R5020-E2_Progesterone, levels = design)
fit2 <- contrasts.fit(fit, contrasts) %>% eBayes()
results_Progesterone <- topTable(fit2, coef = 1, adjust = "BH")
results_R5020 <- topTable(fit2, coef = 2, adjust = "BH")


results <- decideTests(fit2, p.value = 0.01, lfc = 1)
sig <- which(apply(results,1, function(x) sum(x!=0) > 0))


saveRDS(T47D_eset,file = "/gne/research/workspace/schepa/Carroll_Data/T47D_eset.Rds")











