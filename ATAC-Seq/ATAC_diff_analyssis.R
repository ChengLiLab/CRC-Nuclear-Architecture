###perform differential analysis for ATAC-Seq

setwd("/lustre/user/liclab/ganjb/Projects/colon_cancer/ATAC-Seq/analysis_210604/deseq/")
library(data.table)
library(DESeq2)
library(ggpubr)

count <- as.data.frame(fread("./all_peak_expr.txt"))
colnames(count)[1:3] <- c("CHR","START","END")
rownames(count) <- paste0(count$CHR,"_",count$START,"_",count$END)
meta <- as.data.frame(fread("../SampleSheet_all_lihao.csv"))
colnames(count)[4:25] <- meta$SampleID
coldata <- as.data.frame(meta$Treatment)
rownames(coldata) <- meta$SampleID
colnames(coldata) <- "treatment"
head(coldata)
dds <- DESeqDataSetFromMatrix(countData = count[,4:25],colData = coldata,design = ~ treatment)
dds2 <- DESeq(dds)
resultsNames(dds2)
res <-  results(dds2, contrast=c("treatment","normal","tumor"))
res <- as.data.frame(res)
tmp <- res[res$pvalue<=0.05,]
NP_down <- tmp[tmp$log2FoldChange>0,]
NP_down <- data.frame("name"=rownames(NP_down))
NP_down <- tidyr::separate(NP_down,"name",c("CHR","START","END"),"_")
NP_down <- NP_down[is.na(NP_down$END)==F,]
fwrite(NP_down,"NT_down.bed",sep = "\t",col.names = F)
NP_up <- tmp[tmp$log2FoldChange<0,]
NP_up <- data.frame("name"=rownames(NP_up))
NP_up <- tidyr::separate(NP_up,"name",c("CHR","START","END"),"_")
NP_up <- NP_up[is.na(NP_up$END)==F,]
fwrite(NP_up,"NT_up.bed",sep = "\t",col.names = F)
res2 <- results(dds2, contrast=c("treatment","tumor","metastasis"))
res2 <- as.data.frame(res2)
tmp2 <- res2[res2$pvalue<=0.05,]
TM_down <- tmp2[tmp2$log2FoldChange>0,]
TM_down <- data.frame("name"=rownames(TM_down))
TM_down <- tidyr::separate(TM_down,"name",c("CHR","START","END"),"_")
TM_down <- TM_down[is.na(TM_down$END)==F,]
fwrite(TM_down,"TM_down.bed",sep = "\t",col.names = F)
TM_up <- tmp2[tmp2$log2FoldChange<0,]
TM_up <- data.frame("name"=rownames(TM_up))
TM_up <- tidyr::separate(TM_up,"name",c("CHR","START","END"),"_")
TM_up <- TM_up[is.na(TM_up$END)==F,]
fwrite(TM_up,"TM_up.bed",sep = "\t",col.names = F)

###peak expression normalization
rld <- rlogTransformation(dds2) 
exprSet_new=as.data.frame(assay(rld))
par(cex = 0.7)
n.sample=ncol(count[,4:25])
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
#par(mfrow=c(2,2))
boxplot(count[,4:25], col = cols,main="expression value",las=2)
boxplot(exprSet_new, col = cols,main="expression value",las=2)
hist(count[,4:25])
hist(exprSet_new)
options(scipen = 3)
exprSet_new$name <- rownames(exprSet_new)
exprSet_new <- tidyr::separate(exprSet_new,"name",c("CHR","START","END"),"_")
exprSet_new <- exprSet_new[,c(23:25,1:22)]
fwrite(exprSet_new,"All_peak_exp_nor.txt",sep="\t",col.names = F)

###peak num of all sample 
all_num <- as.data.frame(fread("peak_num.txt"))
colnames(all_num) <- c("num","name")
all_num <- tidyr::separate(all_num,"name",c("name","dump"),"_peaks")
barplot(all_num$num,border = NA)
text(cex = 1,x=seq(0.75,1.25*nrow(all_num),1.25),y=-6000,c(all_num$name),xpd = TRUE,srt = 45)
######plot number of diff peaks
data <- data.frame(Tumorigenesis_up=16570,Tumorigenesis_down=14184,
                   Metastasis_up=15719,Metastasis_down=15347)
barplot(as.matrix(data),col="grey",cex.axis=1.3,main="Number of diff peaks",
        border=NA,cex.main=1.5,cex.names=1.3,ylim=c(0,22000),xaxt="n")
text(cex = 1,x=c(0.75,1.95,3.05,4.25),y=-3000,c("T_up","T_down","M_up","M_down"),xpd = TRUE,srt = 45)

#annotate peaks
library(ChIPseeker)
T_up <- readPeakFile("./motif/NT_up_sorted.bed",header=F)
T_down <- readPeakFile("./motif/NT_down_sorted.bed",header=F)
M_up <- readPeakFile("./motif/TM_up_sorted.bed",header=F)
M_down <- readPeakFile("./motif/TM_down_sorted.bed",header=F)
#annotate each cluster
library(GenomicFeatures)
library(clusterProfiler)
library("org.Hs.eg.db")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peaks <- list(Tumorigenesis_up=T_up,Tumorigenesis_down=T_down,Metastasis_up=M_up,Metastasis_down=M_down)
#define promotor region
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter)
#you can transfer geneID (Entrez，ENSEMBL，SYMBOL，GENENAME) after annotatePea be imported into annoDb
peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000), 
                       verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Hs.eg.db")
###visulization
plotAnnoBar(peakAnnoList)
#plotAnnoPie(peakAnnoList$cluster1)
plotDistToTSS(peakAnnoList,title="Distribution of transcription factor-binding loci \n relative to TSS")
