###plot boxplot for single gene in our own dataset and TCGA dataset

setwd("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/onTAD/N4T7_ontad")
library(data.table)
library(tidyr)
library(ggpubr)

func <- function(genename,dir){
  expr <- as.data.frame(fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/RNA-Seq/NvsCTvsLT/20190905/All_read_count_normalized2.txt"))
  apply(expr[,2:10], 2, sum)
  gene_expr <- expr[expr$V1==genename,]
  N <- gene_expr[,c("CRC-01-N","CRC-02-N","CRC-03-N","CRC-04-N","CRC-05-N","CRC-06-N","CRC-07-N","CRC-08-N","CRC-09-N",
                    "CRC-10-N","CRC-11-N","CRC-13-N","CRC-14-N")]
  T <- gene_expr[,c("CRC-01-T","CRC-02-T","CRC-03-T","CRC-04-T","CRC-05-T","CRC-06-T","CRC-07-T","CRC-08-T","CRC-09-T",
                    "CRC-10-T","CRC-11-T","CRC-13-T","CRC-14-T")]
  wilcox.test(as.numeric(N),as.numeric(T))
  #boxplot(list("N"=as.numeric(N),"T"=as.numeric(T)))
  test <- as.data.frame(rbind(t(N),t(T)))
  test$Patient <- "N"
  colnames(test)[1] <- "Expression"
  test[(ncol(N)+1):nrow(test),2] <- "T"
  pdf(paste0(dir,"/",genename,"_own.pdf"),2,3)
  p <- ggboxplot(test, x="Patient",y="Expression",outlier.shape=NA,size=0.5,xlab=F,x.text.angle=0,main=genename,
                 color="Patient",palette=c("lightblue","orange"),add="jitter",add.params=list(size=1.5))
  p <- p + aes(group = Patient)
  p <- p+stat_compare_means()
  print(p)
  dev.off()
  
  ###TCGA data
  ref <- as.data.frame(fread("/lustre/user/liclab/ganjb/resource/hg19_refgene/geneID_ref.txt",select=c(10,2)))
  colnames(ref)<-c("ENSEMBL","SYMBOL")
  expr <- as.data.frame(fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/RNA-Seq/NvsCTvsLT/20190905/TCGA/tcga.txt"))
  for (i in 2:ncol(expr)){
    expr[,i] <- 200000*expr[,i]/sum(expr[,i])
  }
  apply(expr[,2:10],2,sum)
  meta <- as.data.frame(fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/RNA-Seq/NvsCTvsLT/20190905/TCGA/meta_tcga.txt"))
  colnames(expr)[1] <- "ENSEMBL"
  expr$ENSEMBL <- separate(expr,"ENSEMBL",".")[[1]]
  
  clu1_pro <- ref[ref$SYMBOL==genename,]
  clu1_pro <- clu1_pro[is.na(clu1_pro$ENSEMBL)==F,]
  mygene <- as.data.frame(merge(clu1_pro,expr,"ENSEMBL"))
  tcga_N <- mygene[,paste0(c(meta[meta$tissue=="Normal",1]))]
  tcga_T <- mygene[,paste0(c(meta[meta$tissue=="Tumor",1]))]
  wilcox.test(as.numeric(tcga_N),as.numeric(tcga_T))
  #boxplot(list("N"=as.numeric(tcga_N),"T"=as.numeric(tcga_T)))
  test <- as.data.frame(rbind(t(tcga_N),t(tcga_T),use.names=F))
  test$Patient <- "N"
  test[(ncol(tcga_N)+1):nrow(test),2] <- "T"
  colnames(test)[1] <- "Expression"
  pdf(paste0(dir,"/",genename,"_TCGA.pdf"),2,3)
  p <- ggboxplot(test, x="Patient",y="Expression",outlier.shape=NA,size=0.5,xlab=F,x.text.angle=0,main=genename,
                 color="Patient",palette=c("lightblue","orange"),add="jitter",add.params=list(size=0.5))
  p <- p + aes(group = Patient)+stat_compare_means()
  print(p)
  dev.off()
}
dir <- "other"
genename <- "MYCN"
func(genename,dir)

##M-T
genename <- "NOTCH1"
expr <- as.data.frame(fread("E:/ganjingbo/colon-cancer/colon_cancer/RNA-Seq/analysis/All_read_count_normalized2.txt"))
apply(expr[,2:10], 2, sum)
gene_expr <- expr[expr$V1==genename,]
T <- gene_expr[,c("CRC-02-T","CRC-03-T","CRC-04-T","CRC-05-T","CRC-13-T","CRC-14-T")]
M <- gene_expr[,c("CRC-02-M","CRC-03-M","CRC-04-M","CRC-05-M","CRC-13-M","CRC-14-M")]
#boxplot(list("T"=as.numeric(T),"M"=as.numeric(M)))
#T <- as.data.frame(t(T));M<-as.data.frame(t(M))
wilcox.test(as.numeric(T),as.numeric(M),paired = F)

test <- as.data.frame(rbind(t(T),t(M)))
test$Patient <- "T"
colnames(test)[1] <- "Expression"
test[(ncol(T)+1):nrow(test),2] <- "M"
pdf(paste0(genename,"_own_TM.pdf"),2,3)
p <- ggboxplot(test, x="Patient",y="Expression",outlier.shape=NA,size=0.5,xlab=F,x.text.angle=0,main=genename,
               color="Patient",palette=c("orange","pink"),add="jitter",add.params=list(size=1.5))
p + aes(group = Patient)
wilcox.test(as.numeric(M),as.numeric(T))
dev.off()
