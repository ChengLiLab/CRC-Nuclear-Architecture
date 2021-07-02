###analysis the expression changes in AB associated genes

setwd("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/matrix_40k/AB_Compartments")
library(data.table)
library(tidyr)
library(dplyr)
library(ggpubr)
options(scipen = 200)

###gene annotation for AB switch regions
##bedtools sort -i AB_region.bed > AB_region_sorted.bed 
##bedtools intersect -a AB_region_sorted.bed -b /lustre/user/liclab/ganjb/resource/hg19_refgene/TSS_hg19_sorted.txt -wa -wb > AB_region_genes.txt

##TCGA dataset
clu1_pro <- as.data.frame(fread("./diff/test_A2B_NT_genes.txt",select = c(8:9)))
colnames(clu1_pro)<-c("STRAND","SYMBOL")
ref <- as.data.frame(fread("/lustre/user/liclab/ganjb/resource/hg19_refgene/geneID_ref.txt",select=c(10,2)))
colnames(ref)<-c("ENSEMBL","SYMBOL")
clu1_pro <- merge(clu1_pro,ref,"SYMBOL")
clu1_pro <- clu1_pro[,-2]

expr <- as.data.frame(fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/RNA-Seq/NvsCTvsLT/20190905/TCGA/tcga.txt"))
for (i in 2:ncol(expr)){
  expr[,i] <- 200000*expr[,i]/sum(expr[,i])
}
apply(expr[,2:10],2,sum)
meta <- as.data.frame(fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/RNA-Seq/NvsCTvsLT/20190905/TCGA/meta_tcga.txt"))
colnames(expr)[1] <- "ENSEMBL"
expr$ENSEMBL <- separate(expr,"ENSEMBL",".")[[1]]
mygene <- as.data.frame(merge(clu1_pro,expr,"ENSEMBL"))
tcga_N <- mygene[,paste0(c(meta[meta$tissue=="Normal",1]))]
tcga_T <- mygene[,paste0(c(meta[meta$tissue=="Tumor",1]))]
tcga_N$x <- apply(tcga_N,1,mean)
tcga_T$x <- apply(tcga_T,1,mean)
wilcox.test(tcga_N$x,tcga_T$x,paired = T,alternative = "less")
boxplot(list("N"=tcga_N$x,"T"=tcga_T$x))

test <- as.data.frame(rbind(as.matrix(tcga_N[,51:52]),as.matrix(tcga_T[,637:638]),use.names=F))
test$Patient <- "N"
test$Expression <- test$x
test <- test[,-c(1:2)]
test[(nrow(tcga_N)+1):nrow(test),1] <- "T"
p <- ggboxplot(test, x="Patient",y="Expression",outlier.shape=NA,size=0.5,xlab=F,x.text.angle=0,
               color="Patient",palette=c("lightblue","orange"),add="jitter",add.params=list(size=0.5))
p + stat_compare_means(aes(group = Patient),label="p.signif",method="wilcox.test",paired = T)

##our own dataset
clu1_pro <- as.data.frame(fread("./diff/test_B2A_TM_genes.txt",select = c(8:9)))
colnames(clu1_pro)<-c("STRAND","SYMBOL")
ref <- as.data.frame(fread("/lustre/user/liclab/ganjb/resource/hg19_refgene/geneID_ref.txt",select=c(10,2)))
colnames(ref)<-c("ENSEMBL","SYMBOL")
clu1_pro <- merge(clu1_pro,ref,"SYMBOL")
clu1_pro <- clu1_pro[,-2]

expr <- fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/RNA-Seq/NvsCTvsLT/All_read_count_normalized.txt")
apply(expr[,2:10], 2, sum)
gene_expr <- match(clu1_pro$SYMBOL,expr$V1)
gene_expr <- expr[gene_expr,]
gene_expr <- gene_expr[is.na(gene_expr$V1)!=T,]

##tumor vs normal
CT <- gene_expr[,c(15:21,24,27)]
CT$x <- apply(CT,1,mean)
CN <- gene_expr[,c(2:8,10,14)]
CN$x <- apply(CN,1,mean)
wilcox.test(CT$x,CN$x,paired = T)
t.test(CT$x,CN$x,paired = T)
boxplot(list("N"=CN$x,"T"=CT$x))
test <- rbind(CN,CT,use.names=FALSE)
test$Patient <- "N"
test$Expression <- test$x
test <- test[,-c(1:10)]
test[(nrow(CN)+1):nrow(test),1] <- "T"
p <- ggboxplot(test, x="Patient",y="Expression",outlier.shape=NA,size=0.5,xlab=F,x.text.angle=0,
               color="Patient",palette=c("lightblue","orange"),add="jitter",add.params=list(size=0.5))
p + stat_compare_means(aes(group = Patient),label="p.signif",method="wilcox.test",paired = T)

##metastasis vs tumor
CT <- gene_expr[,c(16:19,26,27)]
CT$x <- apply(CT,1,sum);CT$x <- CT$x/6
LT <- gene_expr[,c(29:34)]
LT$x <- apply(LT,1,sum);LT$x <- LT$x/6
wilcox.test(CT$x,LT$x,paired = T,alternative = "greater")
boxplot(list("T"=CT$x,"M"=LT$x))
test <- rbind(CT,LT,use.names=FALSE)
test$Patient <- "T"
test$Expression <- test$x
test <- test[,-c(1:7)]
test[(nrow(CT)+1):nrow(test),1] <- "M"
p <- ggboxplot(test, x="Patient",y="Expression",outlier.shape=NA,size=0.5,xlab=F,x.text.angle=0,
               color="Patient",palette=c("orange","pink"),add="jitter",add.params=list(size=0.5))
p + stat_compare_means(aes(group = Patient),label="p.signif",method="wilcox.test",paired = T)


## GO enrichment by clusterProfiler
gene <- list("B2A_TM"=clu1_pro$ENSEMBL)
library(clusterProfiler)
library("org.Hs.eg.db")
ego <- enrichGO(gene = gene[[1]], OrgDb = org.Hs.eg.db,keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
ego <- as.data.frame(ego@result)
fwrite(ego,"./diff/B2A_TM_genes_GO.csv",sep=",",quote=F)
#dotplot(ego, showCategory=10)
############################## plot GO terms by myself
options(scipen = 200)

##plot top 20 TFs by p-value
tf <- as.data.frame(fread("./diff/B2A_TM_genes_GO.csv",select = c(2,5),nrows = 25))
colnames(tf) <- c("GO Terms","P-value")
tf <- tf[as.numeric(tf$`P-value`)<=0.01,]
tf$`P-value` <- -log10(tf$`P-value`)
tf <- tf[1:20,]
rownames(tf) <- tf$`GO Terms`
tf <- tf[order(as.numeric(tf$`P-value`)),]
barplot(as.numeric(tf$`P-value`),horiz=T,xlim=c(0,1.2*max(tf$`P-value`)),axes=T,col="pink",xlab ="-log10(p-value)",
        cex.axis=1.3,cex.lab=1.5,border = NA) 
for (i in 1:nrow(tf)){
  text(0,(1.2*i-0.6),tf$`GO Terms`[i],cex=1.2,pos=4)
}
