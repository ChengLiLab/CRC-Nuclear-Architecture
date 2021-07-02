###peak-gene association analysis
##step1: annotate peak associated genes (shell): 
##bedtools sort -i peaks.bed > peaks_sorted.bed
##bedtools intersect -a peaks_sorted.bed -b hg19_TSS.txt -wa -wb > peaks_genes.txt

##step2: comparision of peak associated genes expression in different tissues
setwd("/lustre/user/liclab/ganjb/Projects/colon_cancer/ATAC-Seq/analysis_210604/deseq/")
library(data.table)
library(ggpubr)

###gene expression -- TCGA
clu1_pro <- as.data.frame(fread("./NT_down_genes.txt",select = c(8:9)))
colnames(clu1_pro)<-c("STRAND","SYMBOL")
ref <- as.data.frame(fread("/lustre/user/liclab/ganjb/resource/hg19_refgene/geneID_ref.txt",select=c(10,2)))
colnames(ref)<-c("ENSEMBL","SYMBOL")
clu1_pro <- merge(clu1_pro,ref,"SYMBOL")
clu1_pro <- clu1_pro[,-2]

expr <- as.data.frame(fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/RNA-Seq/NvsCTvsLT/20190905/TCGA/tcga_normalized.txt"))
for (i in 2:ncol(expr)){
  expr[,i] <- 200000*expr[,i]/sum(expr[,i])
}
apply(expr[,2:10],2,sum)
meta <- as.data.frame(fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/RNA-Seq/NvsCTvsLT/20190905/TCGA/meta_tcga.txt"))
colnames(expr)[1] <- "ENSEMBL"
expr$ENSEMBL <- tidyr::separate(expr,"ENSEMBL",".")[[1]]
mygene <- as.data.frame(merge(clu1_pro,expr,"ENSEMBL"))
tcga_N <- mygene[,paste0(c(meta[meta$tissue=="Normal",1]))]
tcga_T <- mygene[,paste0(c(meta[meta$tissue=="Tumor",1]))]
tcga_N$x <- apply(tcga_N,1,mean)
tcga_T$x <- apply(tcga_T,1,mean)
wilcox.test(tcga_N$x,tcga_T$x,paired = T,alternative = "less")
boxplot(list("N"=tcga_N$x,"T"=tcga_T$x))

test <- as.data.frame(rbind(as.matrix(tcga_N[,51:52]),as.matrix(tcga_T[,637:638])))
test$Patient <- "N"
test$Expression <- test$x
test <- test[,-c(1:2)]
test[(nrow(tcga_N)+1):nrow(test),1] <- "T"
means <- t(data.frame(mean(tcga_N$x),mean(tcga_T$x)))
colnames(means) <- c("mean")
means<-as.data.frame(means)
means$mean <- round(means$mean,2)
p <- ggboxplot(test, x="Patient",y="Expression",outlier.shape=NA,size=0.5,xlab=F,x.text.angle=0,
               color="Patient",palette=c("lightblue","orange"),add="NA",add.params=list(size=0.5))+
  stat_compare_means(data=test,aes(group = Patient),label="p.format",method="wilcox.test",paired = T,label.x = 1.5)+
  geom_text(data=means, aes(label = mean,x = c(1,2), y = 1.05*mean ))
p

###gene expression -- own
##Tumor vs Normal
clu1_pro <- as.data.frame(fread("./NT_up_genes.txt",select = c(8:9)))
colnames(clu1_pro)<-c("STRAND","SYMBOL")
ref <- as.data.frame(fread("/lustre/user/liclab/ganjb/resource/hg19_refgene/geneID_ref.txt",select=c(10,2)))
colnames(ref)<-c("ENSEMBL","SYMBOL")
clu1_pro <- merge(clu1_pro,ref,"SYMBOL")
clu1_pro <- clu1_pro[,-2]

expr <- fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/RNA-Seq/NvsCTvsLT/All_read_count_normalized.txt")
#expr <- fread("all_sample_gene_expr.txt")
apply(expr[,2:10], 2, sum)
gene_expr <- match(clu1_pro$SYMBOL,expr$V1)
gene_expr <- expr[gene_expr,]
gene_expr <- gene_expr[is.na(gene_expr$V1)!=T,]
CT <- gene_expr[,c(15:21,24,27)]
CT$x <- apply(CT,1,mean)
CN <- gene_expr[,c(2:8,10,14)]
CN$x <- apply(CN,1,mean)
wilcox.test(CT$x,CN$x,paired = T,alternative = "less")
t.test(CT$x,CN$x,paired = T)
boxplot(list("N"=CN$x,"T"=CT$x))
test <- rbind(CN,CT,use.names=FALSE)
test$Patient <- "N"
test$Expression <- test$x
test <- test[,-c(1:10)]
test[(nrow(CN)+1):nrow(test),1] <- "T"
means <- t(data.frame(mean(CN$x),mean(CT$x)))
colnames(means) <- c("mean")
means<-as.data.frame(means)
means$mean <- round(means$mean,2)
p <- ggboxplot(test, x="Patient",y="Expression",outlier.shape=NA,size=0.5,xlab=F,x.text.angle=0,
               color="Patient",palette=c("lightblue","orange"),add="NA",add.params=list(size=0.5))+
  stat_compare_means(data=test,aes(group = Patient),label="p.format",method="wilcox.test",paired = T,label.x = 1.5)+
  geom_text(data=means, aes(label = mean,x = c(1,2), y = 1.05*mean ))
p

##Metastasis vs Tumor
clu1_pro <- as.data.frame(fread("./TM_up_genes.txt",select = c(8:9)))
colnames(clu1_pro)<-c("STRAND","SYMBOL")
ref <- as.data.frame(fread("/lustre/user/liclab/ganjb/resource/hg19_refgene/geneID_ref.txt",select=c(10,2)))
colnames(ref)<-c("ENSEMBL","SYMBOL")
clu1_pro <- merge(clu1_pro,ref,"SYMBOL")
clu1_pro <- clu1_pro[,-2]

expr <- fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/RNA-Seq/NvsCTvsLT/All_read_count_normalized.txt")
#expr <- fread("all_sample_gene_expr.txt")
apply(expr[,2:10], 2, sum)
gene_expr <- match(clu1_pro$SYMBOL,expr$V1)
gene_expr <- expr[gene_expr,]
gene_expr <- gene_expr[is.na(gene_expr$V1)!=T,]
CT <- gene_expr[,c(16:19,26,27)]
CT$x <- apply(CT,1,mean)
LT <- gene_expr[,c(29:34)]
LT$x <- apply(LT,1,mean)
wilcox.test(CT$x,LT$x,paired = T,alternative = "greater")
boxplot(list("T"=CT$x,"M"=LT$x))
test <- rbind(CT,LT,use.names=FALSE)
test$Patient <- "T"
test$Expression <- test$x
test <- test[,-c(1:7)]
test[(nrow(CT)+1):nrow(test),1] <- "M"
means <- t(data.frame(mean(CT$x),mean(LT$x)))
colnames(means) <- c("mean")
means<-as.data.frame(means)
means$mean <- round(means$mean,2)
p <- ggboxplot(test, x="Patient",y="Expression",outlier.shape=NA,size=0.5,xlab=F,x.text.angle=0,
               color="Patient",palette=c("orange","pink"),add="NA",add.params=list(size=0.5))+
  stat_compare_means(data=test,aes(group = Patient),label="p.format",method="wilcox.test",paired = T,label.x = 1.5)+
  geom_text(data=means, aes(label = mean,x = c(1,2), y = 1.05*mean ))
p
