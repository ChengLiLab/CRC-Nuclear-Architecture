## Plot AB compartments after AB_Call.R run
## step1: plot Whole Genome AB compartment
## step2: plot chr1 AB

library(data.table)
setwd("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/merged_matrix_40k/AB_Compartments")

### step1 plot Whole Genome AB compartment
## calculate length for each chromosome
seqname <- fread("CN_merged_40k_compartments.txt",select = 3)
len <- data.frame(chr=character(),len=numeric(),stringsAsFactors = F)
for (i in c(1:22,"X")) {
  chr <- paste0("chr",i)
  pat <- paste0("^",chr,"$") ##used for fully match
  if(i=="X"){i <- 23}
  if(i == 1){
    j = 1
    len[j,1] = chr
    len[j,2] = length(grep(pat,seqname$seqnames))
  }else{
    j = j+1
    len[j,1] <- chr
    len[j,2] <- len[j-1,2]+length(grep(pat,seqname$seqnames))
  }
}

## input data
func.plot <- function(tissue){
  compartments <- fread(paste0(tissue,"_merged_40k_compartments.txt"))
  len <- data.frame(chr=character(),len=integer(),stringsAsFactors=FALSE)
  score<-compartments$score
  score[is.na(score)] <- 0
  return(score)
}
res <- lapply(c("CN","CT","LT"), func.plot)

## start plotting
png("AB_WG.png",width = 1000,height = 520)
par(mfrow=c(3,1),mar=c(0,0,0,8),oma=c(4,4,3,0),mgp=c(3,2,0),xpd=T)
barplot(res[[1]],col = ifelse(res[[1]]>0, "darkred", "darkblue"),border = NA,axes = F,ylim = c(-1.1*max(res[[1]]),1.1*max(res[[1]])),
        space=0,cex.axis = 2)
mtext("CN",side = 2,cex = 2,las=2)
mtext("Whole Genome AB Compartments Distribution",side = 3,cex=2.5)
barplot(res[[2]],col = ifelse(res[[2]]>0, "darkred", "darkblue"),border = NA,axes = F,ylim = c(-1.1*max(res[[1]]),1.1*max(res[[1]])),
        space=0,cex.axis = 2)
mtext("CT",side = 2,cex = 2,las=2)
legend(length(res[[1]]),max(res[[1]]),c("A","B"),col=c("darkred","darkblue"),bty="n", pch=15,cex = 4,xpd=T)
barplot(res[[3]],col = ifelse(res[[3]]>0, "darkred", "darkblue"),border = NA,axes = F,ylim = c(-1.1*max(res[[1]]),1.1*max(res[[1]])),
        space=0,cex.axis = 2)
mtext("LT",side = 2,cex = 2,las=2)
lines(0,0,lwd=3)
axis(1, at=c(0,len$len),labels=c("","chr1","","chr3","","chr5","","chr7","","chr9","","chr11","","chr13","","chr15","","chr17",
                                 "","chr19","","chr21","","chrX"),cex.axis = 2.5)
dev.off()


### step2: plot chr1 AB
func.plot <- function(tissue){
  compartments <- fread(paste0(tissue,"_merged_40k_compartments.txt"))
  compartments <- compartments[grep("^chr1$",compartments$seqnames)]
  score<-compartments$score
  score[is.na(score)] <- 0
  return(score)
}
res <- lapply(c("CN","CT","LT"), func.plot)

png("AB_chr1.png",width = 1000,height = 520)
par(mfrow=c(3,1),mar=c(2,2,0,8),oma=c(4,4,5,0),mgp=c(3,2,0),xpd=T)
barplot(res[[1]],col = ifelse(res[[1]]>0, "darkred", "darkblue"),border = NA,axes = F,ylim = c(-1.1*max(res[[1]]),1.1*max(res[[1]])),
        space=0,cex.axis = 2)
mtext("CN",side = 2,cex = 2,las=2)
mtext("Chromosome 1 AB Compartments Distribution",side = 3,cex=2.5)
barplot(res[[2]],col = ifelse(res[[2]]>0, "darkred", "darkblue"),border = NA,axes = F,ylim = c(-1.1*max(res[[1]]),1.1*max(res[[1]])),
        space=0,cex.axis = 2)
mtext("CT",side = 2,cex = 2,las=2)
legend(length(res[[1]]),max(res[[1]]),c("A","B"),col=c("darkred","darkblue"),bty="n", pch=15,cex = 4,xpd=T)
barplot(res[[3]],col = ifelse(res[[3]]>0, "darkred", "darkblue"),border = NA,axes = F,ylim = c(-1.1*max(res[[1]]),1.1*max(res[[1]])),
        space=0,cex.axis = 2)
mtext("LT",side = 2,cex = 2,las=2)
axis(1, at=c(seq(0,6000,1000),length(res[[1]])),labels = c("0","100M","200M","300M","400M","500M","600M",""),lwd = 3,cex.axis=2)
dev.off()
