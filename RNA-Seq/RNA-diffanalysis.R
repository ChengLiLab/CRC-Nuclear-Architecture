###perform paired differential analysis for RNA-Seq data and plot volcano plot 

setwd("/lustre/user/liclab/ganjb/Projects/colon_cancer/RNA-Seq/NvsCTvsLT/20190905")
library(data.table)

####Step1: normalize expression data and do diff express analysis by deseq2
data <- data.frame(fread("./All_read_count_plus_public.txt"))
rownames(data) <- data$Gene_Symbol
data <- data[,-1]

library(DESeq2)
colData <- as.matrix(fread("../rna-metadata.CSV",header = T,sep = "," ,select=2))
name_colData <- fread("../rna-metadata.CSV",header = T,sep = "," ,select=1)
rownames(colData) <- as.matrix(name_colData)
colData <- as.data.frame(colData)
head(colData)
data <- as.data.frame(data)
dds <- DESeqDataSetFromMatrix(data, colData, design= ~ treatment)
###normlization
# rld <- rlogTransformation(dds) #slow than vst()
rld <- vst(dds)
data_new=assay(rld)
#par(cex = 0.7)
n.sample=ncol(data)
#if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
par(mar=c(7,5,2,2))
boxplot(data, col = cols,main="",ylab="Raw Expression",las=2,cex.lab =1.5,outline=F)
boxplot(data_new, col = cols,main="",ylab="Normalized Expression",las=2,cex.lab =1.5)
write.table(data,file = "All_read_count_new.txt",quote = F)
write.table(data_new,file = "All_read_count_normalized_new.txt",quote = F)
###diff analysis
dds2 <- DESeq(dds)
resultsNames(dds2)

res1 <-  results(dds2, contrast=c("treatment","pri_tumor","normal"))
PvsN <- as.data.frame(res1)
fwrite(PvsN,"NP.txt",sep="\t",quote=F,row.names = T)
PvsN_down <- down <- PvsN[PvsN$pvalue<=0.05 & PvsN$log2FoldChange<=-1 & is.na(PvsN$pvalue)==F,]
fwrite(as.data.frame(row.names(down)),"NP_down.name")
PvsN_up <- up <- PvsN[PvsN$pvalue<=0.05 & PvsN$log2FoldChange>=1 & is.na(PvsN$pvalue)==F,]
fwrite(as.data.frame(row.names(up)),"NP_up.name")

res2 <-  results(dds2, contrast=c("treatment","liver_tumor","pri_tumor"))
MvsP <- as.data.frame(res2)
fwrite(MvsP,"PM.txt",sep="\t",quote=F,row.names=T)
MvsP_down <- down <- MvsP[MvsP$pvalue<=0.05 & MvsP$log2FoldChange<=-1 & is.na(MvsP$pvalue)==F,]
fwrite(as.data.frame(row.names(down)),"PM_down.name")
MvsP_up <- up <- MvsP[MvsP$pvalue<=0.05 & MvsP$log2FoldChange>=1 & is.na(MvsP$pvalue)==F,]
fwrite(as.data.frame(row.names(up)),"PM_up.name")

res3 <-  results(dds2, contrast=c("treatment","liver_tumor","normal"))
MvsN <- as.data.frame(res3)
fwrite(MvsN,"NM.txt",sep="\t",quote=F,row.names=T)
MvsN_down <- down <- MvsN[MvsN$pvalue<=0.05 & MvsN$log2FoldChange<=-1 & is.na(MvsN$pvalue)==F,]
fwrite(as.data.frame(row.names(down)),"NM_down.name")
MvsN_up <- up <- MvsN[MvsN$pvalue<=0.05 & MvsN$log2FoldChange>=1 & is.na(MvsN$pvalue)==F,]
fwrite(as.data.frame(row.names(up)),"NM_up.name")

NPM <- unique(c(row.names(PvsN_down),row.names(MvsP_down),row.names(MvsN_down),row.names(PvsN_up),
              row.names(MvsP_up),row.names(MvsN_up)))
library(VennDiagram)
venn.diagram(list(T=row.names(PvsN_down),M=row.names(MvsN_down)), fill=c("orange","pink"), alpha=c(0.5,0.5),
             imagetype="png",filename="Venn_down.png",col="transparent",cex=3,cat.cex=3)
venn.diagram(list(P=row.names(PvsN_up),M=row.names(MvsN_up)), fill=c("orange","pink"), alpha=c(0.5,0.5), 
             imagetype="png",filename="Venn_up.png",col="transparent",cex=3,cat.cex=3)
intersect(row.names(PvsN_down),row.names(MvsN_down))
length(intersect(row.names(PvsN_down),row.names(MvsN_up)))

#### Step2: cluster all samples and plot heatmap of top variant genes expression
###clustering
#sd <- apply(data_new,1,sd)
#sd_1000 <- data_new[order(sd,decreasing = T)[1:1000],]
#hc <- hclust(dist(t(log(sd_1000+1))))
#plot(hc, hang = -1, cex=1.5, main="",cex.main = 1.5, cex.lab = 1.5)

###heatmap of most variant genes expression
#library(pheatmap)
#annotation_col <- data.frame(rep(c("Normal","Tumor","Metastasis"),times=c(13,13,6)))
#colnames(annotation_col)<-c("Tissues")
#ann_color <- list(Tissues=c(Normal="lightblue",Tumor="Orange",Metastasis="pink"))
#rownames(annotation_col) <- colnames(data_new)
#pheatmap(sd_1000,border_color = NA,show_colnames = T,show_rownames = F, cluster_cols = T,main = "",
#       fontsize = 15, annotation_col = annotation_col,annotation_colors = ann_color,annotation_names_col = F)

#### Step3: volcano plot, MA plot 
###volcano plot NP
file <- "NP.txt"
all <- data.frame(fread(file,sep="\t",select=c(1,3,6,7)))
all <- all[is.na(all$pvalue)==F,]
rownames(all) <- all[,1]
all <- all[,-1]
up <- all[all$log2FoldChange>=log2(1.5) & all$pvalue<=0.05,]
down <- all[all$log2FoldChange<=log2(1/1.5) & all$pvalue<=0.05,]
plot(all$log2FoldChange,-log10(all$pvalue),pch=20,col="black",lwd=1,ylim=c(0,25),
     xlab="log2(T/N)",ylab="-log10(p-value)",cex.lab =1.5,asp=0,cex.axis=1.3,bty="l")
points(down$log2FoldChange,-log10(down$pvalue),pch=20,col="steelblue",asp=0,cex=1)
points(up$log2FoldChange,-log10(up$pvalue),pch=20,col="hotpink",asp=0,cex=1)
text(-6,25,nrow(down),cex=1.3,col="steelblue")
text(6,25,nrow(up),cex=1.3,col="hotpink")
gene <- "MUC2"
x <- all[gene,1]
y <- -log10(all[gene,2])
points(x,y,pch=20,col="red",asp=0,cex=1)
lines(c(x,x),c(y,y+13),type = "l")
text(x,y+14,gene,cex=1.3)
gene <- "SMAD4"
x <- all[gene,1]
y <- -log10(all[gene,2])
points(x,y,pch=20,col="red",asp=0,cex=1)
lines(c(x,x),c(y,y+13),type = "l")
text(x,y+14,gene,cex=1.3)
gene <- "SPP1"
x <- all[gene,1]
y <- -log10(all[gene,2])
points(x,y,pch=20,col="red",asp=0,cex=1)
lines(c(x,x+0.5),c(y,y+13),type = "l")
text(x+0.5,y+14,gene,cex=1.3)
gene <- "MYC"
x <- all[gene,1]
y <- -log10(all[gene,2])
points(x,y,pch=20,col="red",asp=0,cex=1)
lines(c(x,x-1),c(y,y+10),type = "l")
text(x-1,y+11,gene,cex=1.3)
gene <- "FOXP3"
x <- all[gene,1]
y <- -log10(all[gene,2])
points(x,y,pch=20,col="red",asp=0,cex=1)
lines(c(x,x-1),c(y,y+16),type = "l")
text(x-1,y+17,gene,cex=1.3)

###volcano plot MN
file <- "PM.txt"
all <- data.frame(fread(file,sep="\t",select=c(1,3,6,7)))
all <- all[is.na(all$pvalue)==F,]
rownames(all) <- all[,1]
all <- all[,-1]
up <- all[all$log2FoldChange>=log2(1.5) & all$pvalue<=0.05,]
down <- all[all$log2FoldChange<=log2(1/1.5) & all$pvalue<=0.05,]
plot(all$log2FoldChange,-log10(all$pvalue),pch=20,col="black",lwd=1,
     xlab="log2(M/T)",ylab="-log10(p-value)",cex.lab =1.5,asp=0,cex.axis=1.3,bty="l")
points(down$log2FoldChange,-log10(down$pvalue),pch=20,col="steelblue",asp=0,cex=1)
points(up$log2FoldChange,-log10(up$pvalue),pch=20,col="hotpink",asp=0,cex=1)
text(-6,25,nrow(down),cex=1.3,col="steelblue")
text(6,25,nrow(up),cex=1.3,col="hotpink")
gene <- "BMP3"
x <- all[gene,1]
y <- -log10(all[gene,2])
points(x,y,pch=20,col="red",asp=0,cex=1)
lines(c(x,x),c(y,y+13),type = "l")
text(x,y+14,gene,cex=1.3)
gene <- "CASP3"
x <- all[gene,1]
y <- -log10(all[gene,2])
points(x,y,pch=20,col="red",asp=0,cex=1)
lines(c(x,x),c(y,y+13),type = "l")
text(x,y+14,gene,cex=1.3)
gene <- "SPP1"
x <- all[gene,1]
y <- -log10(all[gene,2])
points(x,y,pch=20,col="red",asp=0,cex=1)
lines(c(x,x+1.5),c(y,y+9),type = "l")
text(x+1.5,y+10,gene,cex=1.3)
gene <- "MUC6"
x <- all[gene,1]
y <- -log10(all[gene,2])
points(x,y,pch=20,col="red",asp=0,cex=1)
lines(c(x,x-1.5),c(y,y+8),type = "l")
text(x-1.5,y+9,gene,cex=1.3)
gene <- "TCF4"
x <- all[gene,1]
y <- -log10(all[gene,2])
points(x,y,pch=20,col="red",asp=0,cex=1)
lines(c(x,x),c(y,y+9),type = "l")
text(x,y+10,gene,cex=1.3)

###MA plot
#plotMA(res1)
###PCA plot
#library(ggplot2)
#rld <- vst(dds, blind=FALSE)
#pcaData <- plotPCA(rld, intgroup=c("treatment"), returnData=TRUE) ##1000 most variant genes
#percentVar <- round(100 * attr(pcaData, "percentVar"))
#ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=treatment)) +
#  geom_point(size=3) +
#  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#  coord_fixed()
