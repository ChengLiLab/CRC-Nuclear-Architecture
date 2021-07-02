#####cluster our RNA-Seq data with TCGA COAD and READ dataset

setwd("/lustre/user/liclab/ganjb/Projects/colon_cancer/RNA-Seq/NvsCTvsLT/20190905/TCGA")
library(data.table)
library(dplyr)
library(tidyr)
library(ggbiplot)
library(DESeq2)

tcga <- as.data.frame(fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/RNA-Seq/public/tcga/TCGA-COAD.htseq_counts.tsv"))
meta <- data.frame("sample"=colnames(tcga)[2:513])
meta <- separate(meta,"sample",c("id","tissue","batch","tag"),"-")
meta$id <- colnames(tcga)[2:513]
table(meta$tag)
meta[meta$tag=="01A",2]<-meta[meta$tag=="01B",2]<-meta[meta$tag=="01C",2]<-meta[meta$tag=="02A",2]<-"Tumor"
meta[meta$tag=="11A",2] <- "Normal"
meta[meta$tag=="06A",2] <- "Metastasis"
meta$batch <- "COAD"
fwrite(meta,"meta_tcga_COAD.txt",sep="\t",quote=F)

###PCA for single TCGA dataset
pca <- prcomp(t(tcga[1:20000,2:513]),scale.=F)
#ggscreeplot(pca)
ggbiplot(pca,obs.scale=1,var.scale=1,groups=as.factor(meta$tissue[]),ellipse=T,circle=F,
         var.axes=F)+scale_color_discrete(name = '')+
  theme(legend.direction='horizontal',legend.position='top')


###PCA for merged TCGA dataset
coad <- as.data.frame(fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/RNA-Seq/public/tcga/TCGA-COAD.htseq_counts.tsv"))
read <- as.data.frame(fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/RNA-Seq/public/tcga/TCGA-READ.htseq_counts.tsv"))
tcga <- merge(coad,read,"Ensembl_ID")
tcga <- tcga[-c(1:5),]
apply(tcga[,2:30],2,sum)
#fwrite(tcga,"tcga.txt",sep="\t",quote=F)
library(DESeq2)
rownames(tcga) <- tcga$Ensembl_ID
tcga <- tcga[,-1]
meta <- as.data.frame(fread("meta_tcga.txt"))
rownames(meta) <- meta$id
meta <- meta[,-c(1,3)]
head(meta)
tcga <- round(10*tcga)
dds <- DESeqDataSetFromMatrix(tcga, meta, design= ~ tissue)
####### normlization for TCGA dataset
# rld <- rlogTransformation(dds) #slow than vst()
rld <- vst(dds)
tcga_new=assay(rld)
apply(tcga_new[,2:30],2,sum)+
fwrite(tcga,file = "tcga.txt",quote = F)
fwrite(as.data.frame(tcga_new),"tcga_normalized.txt",quote = F,col.names=T,row.names=T,sep="\t")
pca <- prcomp(t(tcga[1:20000,1:689]),scale.=F)
#ggscreeplot(pca)
ggbiplot(pca,obs.scale=1,var.scale=1,groups=as.factor(meta$tissue[]),ellipse=T,circle=F,
         var.axes=F)+scale_color_discrete(name = '')+
  theme(legend.direction='horizontal',legend.position='top')
#optimized 2D PCA plot
data <- data.frame(pca$x)
data <- data[,c(1,2)]
data$tissue <- meta$tissue
data$batch <- meta$batch
p <- ggplot(data, aes(PC1, PC2, shape = batch))
p + geom_point(aes(colour = tissue), size = 2) +
  labs(x = "PC1 (33.7%)", y = "PC2 (8.2%)")+
  theme_classic()

###PCA for our own data
own <- as.data.frame(fread("./own_hg38.txt"))
rownames(own) <- own$Name
own <- own[,-c(1,4,10)]
meta_own <- fread("../../rna-metadata-batch.CSV")
meta_own <- meta_own[19:50,]
own <- own[,c(meta_own$id)]
apply(own[,],2,sum)
#own_n <- as.data.frame(t(own))
#own_n$sum <- apply(own_n,1,sum)
##for (i in 1:nrow(own_n)){
#  own_n[i,] <- own_n[i,]*120000/own_n[i,60484]
#}
#own_n <- as.data.frame(t(own_n))
#apply(own_n,2,sum)
#own_n <- own_n[-60484,]
#own_n$sd <- apply(own_n,1,sd)
#own_n <- own_n[sort(own_n$sd,decreasing=T),-33]
#pca <- prcomp(t(own_n[1:1000,]),scale.=F)
#ggscreeplot(pca)
#ggbiplot(pca,obs.scale=1,var.scale=1,groups=as.factor(meta_own$Tissue[]),ellipse=T,circle=F,
#         var.axes=F)+scale_color_discrete(name = '')+
#  theme(legend.direction='horizontal',legend.position='top')

dds_own <- DESeqDataSetFromMatrix(own, meta_own, design= ~ Tissue)
plotPCA(DESeqTransform(dds_own),intgroup="Tissue",ntop=500)
#vst
rld_own <- rlogTransformation(dds_own) #slow than vst()
#assay(rld_own) <- limma::removeBatchEffect(assay(rld_own), meta_own$batch)
plotPCA(rld_own,intgroup="Tissue",ntop=60483)
vsd_own <- vst(dds_own)
#assay(vsd_own) <- limma::removeBatchEffect(assay(vsd_own), meta_own$batch)
plotPCA(vsd_own,intgroup="Tissue",ntop=60483)
data_new=as.data.frame(assay(rld))
data_new$sd <- apply(data_new,1,sd)
data_new <- data_new[order(data_new$sd,decreasing = T),-33]
#par(cex = 0.7)
n.sample=ncol(data_new)
#if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
par(mar=c(7,5,2,2))
boxplot(own, col = cols,main="",ylab="Raw Expression",las=2,cex.lab =1.5,outline=F)
boxplot(data_new, col = cols,main="",ylab="Normalized Expression",las=2,cex.lab =1.5)
#write.table(own,file = "All_read_count_hg38.txt",quote = F)
write.table(data_new,file = "own_hg38_rld.txt",quote = F,row.names = T)

###PCA for all data (our own data + TCGA dataset)
own <- as.data.frame(fread("own_hg38.txt"))
tcga <- as.data.frame(fread("tcga.txt"))
tcga[,2:690] <- (2^tcga[,2:690])-1
colnames(own)[1] <- "Ensembl_ID"
own <- own[,-c(4,10)]
meta <- as.data.frame(fread("meta_all.txt"))
own <- own[,c("Ensembl_ID",c(meta$id[690:721]))]
all <- merge(tcga, own,"Ensembl_ID")
all[,2:722] <- round(all[,2:722])
dds <- DESeqDataSetFromMatrix(all[,2:722], meta, design= ~ tissue)
plotPCA(DESeqTransform(dds),intgroup="tissue",ntop=50000)
#rld or vst
#rld <- rlogTransformation(dds) #slow than vst()
#assay(rld) <- limma::removeBatchEffect(assay(rld), meta$batch)
#plotPCA(rld,intgroup="Tissue",ntop=1000)
vsd <- vst(dds)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), meta$batch) #remove batch effect for the two datasets
plotPCA(vsd,intgroup="tissue",ntop=1000)
data_new=as.data.frame(assay(vsd))
data_new$sd <- apply(data_new,1,sd)
data_new <- data_new[order(data_new$sd,decreasing = T),-722]
pca <- prcomp(t(data_new[1:1000,1:721]),scale.=F)
#ggscreeplot(pca)
#ggscreeplot(pca)
ggbiplot(pca,obs.scale=1,var.scale=1,groups=as.factor(meta$tissue[]),ellipse=T,circle=F,
         var.axes=F)+scale_color_discrete(name = '')+
  theme(legend.direction='horizontal',legend.position='top')
#optimized 2D PCA plot
data <- data.frame(pca$x)
data <- data[,c(1,2)]
data$tissue <- meta$tissue
data$batch <- meta$batch
p <- ggplot(data, aes(PC1, PC2, shape = batch))
p + geom_point(aes(colour = tissue), size = 2) +
  labs(x = "PC1 (12.9%)", y = "PC2 (9.6%)")+
  theme_classic()

###cluster by most variant genes in our own dataset
own <- as.data.frame(fread("./own_hg38_rld.txt"))
rownames(own) <- own$V1
own <- own[,-c(1)]
meta_own <- fread("../../rna-metadata-batch.CSV")
meta_own <- meta_own[19:50,]
own <- own[,c(meta_own$id)]
apply(own[,],2,sum)
own$sd <- apply(own,1,sd)
own <- own[order(own$sd,decreasing = T),-33]
pca <- prcomp(t(own[1:500,1:32]),scale.=F)
ggscreeplot(pca)
ggbiplot(pca,obs.scale=1,var.scale=1,groups=as.factor(meta_own$Tissue[]),ellipse=T,circle=F,
         var.axes=F)+scale_color_discrete(name = '')+
  theme(legend.direction='horizontal',legend.position='top')
#pca
data_new=as.data.frame(assay(vsd))
rownames(data_new) <- all$Ensembl_ID
pca <- prcomp(t(data_new[rownames(own[1:500,]),1:721]),scale.=F)
ggbiplot(pca,obs.scale=1,var.scale=1,groups=as.factor(meta$tissue[]),ellipse=T,circle=F,
         var.axes=F)+scale_color_discrete(name = '')+
  theme(legend.direction='horizontal',legend.position='top')
#optimized 2D PCA plot
data <- data.frame(pca$x)
data <- data[,c(1,2)]
data$tissue <- meta$tissue
data$batch <- meta$batch
p <- ggplot(data, aes(PC1, PC2, shape = batch))
p + geom_point(aes(colour = tissue), size = 1.5) +
  scale_colour_manual(values=list("Normal"="lightblue","Tumor"="orange","Metastasis"="pink"))+
  scale_shape_manual(values=list("own"=19,"COAD"=15,"READ"=17))+
  labs(x="PC1 (20.4% Variance)",y="PC2 (13.3% Variance)")+
  theme_classic()
###3D PCA plot
library(scatterplot3d)
data_3d <- as.data.frame(pca$x)
meta$shape <- 19
meta[meta$batch=="own",5] <- 19
meta[meta$batch=="COAD",5] <- 15
meta[meta$shape=="READ",5] <- 17
meta$color <- "pink"
meta[meta$tissue=="Normal",6]<-"lightblue"
meta[meta$tissue=="Tumor",6]<-"orange"
with(data_3d, {
  s3d <- scatterplot3d(PC1, PC2, PC3,   
                       pch=c(meta$shape),        # circle color indicates no. of cylinders
                       color=c(meta$color),
                       scale.y=.75,                # scale y axis (reduce by 25%)
                       main="",
                       xlab="PC1 (15.5%)",
                       ylab="PC2 (8.32%)",
                       zlab="PC3 (7.06%)")
  #  s3d.coords <- s3d$xyz.convert(PC1, PC2, PC3)
  #  text(s3d.coords$x, s3d.coords$y,    # x and y coordinates
  #       labels=class$Tissue,      # text to plot
  #       pos=4, cex=.5)                  # shrink text 50% and place to right of points)
})
