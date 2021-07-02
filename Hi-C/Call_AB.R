## call AB compartments for matrix from HiC-Pro
## step1: construct HiTC list object from whole genome contact matrix in triple list format
## step2: call A/B by pca.hic

## step1 construct HiTC object 
library(data.table)
library(HiTC)
#library(refGenome)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
setwd("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/matrix_40k/AB_Compartments")

func <- function(tissue){
  library(data.table)
  library(HiTC)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  matrix <- paste0("../HIC",tissue,"_40000_iced.matrix")
  #matrix <- "../matrix_40k/HIC04_CN_40000_iced.matrix"
  xgi <- ygi <- "./40k_abs.bed"
  hic <- importC(matrix, xgi, ygi,allPairwise=T)
  hiC14 <- extractRegion(hic$chr12chr12, chr="chr12", from=82840000, to=94520000)
  png(paste0(tissue,"_chr12-82840000-94520000_40k_heatmap.png"),width = 1400,height = 200)
  mapC(hiC14, trim.range=.95)
  #mapC(forceSymmetric(HTClist(hiC14)), trim.range=.95)
  dev.off()
  
  ## step2 extract individual chromosomes interaction maps and call A/B comparments
  human.genes = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

  hic.intra <- hic[isIntraChrom(hic)] ## get all cis maps
  hic.intra
  hic.intra <- hic.intra[1:24] # contain chromosome 1-22 and X
  
  pca.intra <- lapply(hic.intra, function(xx){
    pca <- pca.hic(xx, normPerExpected = T, npc = 1, gene.gr = human.genes)
    return(pca)
  })
  
  compartments <- data.frame(group=character(), group_name=character(), seqnames=character(), start=integer(), end=integer(), width=integer(), strand=character(),
                             score=double(), genedens=integer(), ccompartments=character(),stringsAsFactors=FALSE)
  
  for(i in c(1:22,"X")){
    chrom <- paste0("chr",i,"chr",i)
    pca.chrom <- data.frame(pca.intra[chrom])
    pca.chrom <- pca.chrom[which(pca.chrom[,1] == "1"),]
    names(pca.chrom) <- c("group","group_name","seqnames","start","end","width","strand","score","genedens","ccompartments")
    compartments <- rbind(compartments, pca.chrom)
  }
  ## step4 save results
  fwrite(compartments, paste0(tissue,"_merged_40k_compartments.txt"),
              sep = '\t', row.names = F, col.names = T, quote = F)
  
}
cl <- makeCluster(9)
samples <- c("01_CN","01_CT","02_CN","02_CT","02_LT","02_LN","03_CN","03_CT","03_LT","03_NT","04_CN","04_CT","04_LT","05_CN","05_CT","05_LT",
             "06_CN","06_CT","07_CN","07_CT","08_CN","08_CT","09_CN","09_CT","10_CN","10_CT","11_CN","11_CT","13_CN","13_CT",
             "13_LT","14_CN","14_CT","14_LT")
parLapply(cl,samples,func)
stopCluster(cl)

## plot heatmap of chr1
#hiC14 <- extractRegion(hic$chr1chr1, chr="chr1", from=1.8e+07, to=106368584)
#hiC14 <- extractRegion(hic$chr1chr1, chr="chr1", from=1, to=249250621)
#par(mar=c(4,2,4,2))
#mapC(forceSymmetric(HTClist(hiC14)), trim.range=.95,col.pos=c("white", "orange", "red", "black"))
#mapC(forceSymmetric(HTClist(hiC14)), trim.range=.95)
## Data Normalization by Expected number of Counts
#hiC14norm <- normPerExpected(hiC14, method="mean")
#mapC(forceSymmetric(HTClist(hiC14norm)), log.data=TRUE)
## Correlation Map of Chromosome 14
#pm <- getPearsonMap(hiC14)
#mapC(forceSymmetric(HTClist(pm)), maxrange=1, col.pos=c("black","red"), col.neg=c("black","blue"))
## calculate TAD
#di<-directionalityIndex(hiC14)
#barplot(di, col=ifelse(di>0,"darkred","darkgreen"), space = 0,border = NA)
