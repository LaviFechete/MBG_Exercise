.libPaths(c("/content/miniconda/envs/myenv/lib/R/library/", .libPaths())) 

suppressMessages(library(ggplot2))
suppressMessages(library(ChIPseeker))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

message("Loading data...")
TBP <- read.delim("/content/macs3_results/TBP_summits.bed", header=F)
TBP$V1 <- "chr1"
NF_YA <- read.delim("/content/macs3_results/NF-YA_summits.bed", header=F)
NF_YA$V1 <- "chr1"
write.table(TBP, "/content/macs3_results/TBP_renamed.bed", sep="\t", row.names=F, col.names=F, quote=F)
write.table(NF_YA, "/content/macs3_results/NF-YA_renamed.bed", sep="\t", row.names=F, col.names=F, quote=F)
samplefiles <- list("/content/macs3_results/TBP_renamed.bed", "/content/macs3_results/NF-YA_renamed.bed")
names(samplefiles) <- c("TBP", "NF_YA")

message("Plot peaks...")
p1<- covplot(samplefiles$TBP, weightCol="V5", title = "ChIP Peaks over region - TBP")
suppressMessages(ggsave("covplot_TBP.png", plot = p1, width=20, height=5, dpi=300))

p2 <- covplot(samplefiles$NF_YA, weightCol="V5", title = "ChIP Peaks over region - NF-YA")
suppressMessages(ggsave("covplot_NFYA.png", plot = p2, width=20, height=5, dpi=300))

message("Plot Venn diagram...")
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-3000, 3000), verbose=FALSE)

p3 <-plotAnnoBar(peakAnnoList)
suppressMessages(ggsave("plotAnnoBar.png", plot = p3))

genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
png("vennplot.png")
vennplot(genes)
trash <- dev.off()

message("Done!")


