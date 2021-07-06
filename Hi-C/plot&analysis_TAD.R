### Used to plot TAD boundaries and analysis TAD boundaries
### newer version
library(data.table)
library(parallel)
library(plyr)
options(scipen = 200)
setwd("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/matrix_40k/TAD")

ref <- as.data.frame(fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/matrix_40k/40k_abs.bed"))
colnames(ref) <- c("chr","start","end","bin")
ref <- cbind(ref,data.frame(ref=paste0(ref$chr,":",ref$start+1,"-",ref$end+1)))
path <- "/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/matrix_40k/TAD/filtered_TAD/"
file <- sort(list.files(path,pattern = "boundaries.bed"))
boundary <- data.frame(rep(0,nrow(ref)))
for (i in c(1:length(file))){
  name <- gsub("_filtered_boundaries.bed","",file[i])
  mid <- fread(paste0(path,file[i]),select = 2)
  fwrite(mid,"mid.txt",col.names = F)
  mid <- fread("mid.txt",sep = "|",select = 3,header = F)
  mid2 <- data.frame(mid2=rep(0,nrow(ref)))
  colnames(mid2) <- name
  mid2[match(mid$V3,ref$ref),1] <- 1
  boundary <- cbind(boundary,mid2)
  rm(mid);rm(mid2)
}
boundary <- boundary[,-1]
boundary <- cbind(ref[,1:4],boundary)
fwrite(boundary,"boundary_all_samples.bed",sep="\t")

### start plotting
boundary <- fread("boundary_all_samples.bed",quote=F)
chr7 <- as.data.frame(boundary[boundary$chr=="chr7",])

png("TAD_boundaries_chr7_new.png",width = 1000,height = 500)
par(mfrow=c(25,1),mar=c(0,4.5,0.2,0),oma=c(0,0,0,0),mgp=c(1,0,0),xpd=NA)
for (i in c(5:ncol(chr7))) {
  #max <- round(max(abs(mat_chr2[,i])),1)
  barplot(chr7[,i],col = "darkblue",border = ifelse(chr7[,i]==1,"darkblue","white"),axes = F,ylim = c(0,1),space=0,cex.axis = 2)
  mtext(names(chr7)[i],side = 2,cex = 1,las=2,line = -3)
}
dev.off()

### start analysis
boundary$start <- boundary$start+1
del <- c()
k=1
for (i in c(1:nrow(boundary))){
  test <- as.data.frame(table(boundary[i,] == 0))
  condi <- test[test$Var1=="FALSE",2]
  if(condi == 4){
    del[k] <- i
    k=k+1
  }
}
boundary <- boundary[-del,]

## merge close bin
fwrite(as.data.frame(boundary$chr),"mid.txt",col.names = F)
chr <- fread("mid.txt",sep = "r",header = F,select = 2)
colnames(chr) <- "chr"
chr[chr=="X",1] <- "23"
boundary$chr <- as.data.frame(chr)
boundary$chr <- as.numeric(chr$chr)
del <- boundary[1,]
del <- del[-1,]
i=2
while (i <= nrow(boundary)) {
  if((boundary[i,4]-boundary[i-1,4])<=6){
    del <- rbind(del,ceiling(boundary[i,]+boundary[i-1,])/2);i=i+2
  }else {del <- rbind(del,boundary[i-1,])
  i=i+1
  }
}
del$bin <- ceiling(del$bin)
del[,5:29] <- ceiling(del[,5:29])

del2 <- del[1,]
del2 <- del2[-1,]
i=2
while (i <= nrow(del)) {
  if((del[i,4]-del[i-1,4])<=6){
    del2 <- rbind(del2,ceiling(del[i,]+del[i-1,])/2);i=i+2
  }else {del2 <- rbind(del2,del[i-1,]);i=i+1}
}
del2$bin <- ceiling(del2$bin)
del2[,5:29] <- ceiling(del2[,5:29])


del3 <- del2[1,]
del3 <- del3[-1,]
i=2
while (i <= nrow(del2)) {
  if((del2[i,4]-del2[i-1,4])<=6){
    del3 <- rbind(del3,ceiling(del2[i,]+del2[i-1,])/2);i=i+2
  }else {del3 <- rbind(del3,del2[i-1,]);i=i+1}
}
del3$bin <- ceiling(del3$bin)
del3[,5:29] <- ceiling(del3[,5:29])


del4 <- del3[1,]
del4 <- del4[-1,]
i=2
while (i <= nrow(del3)) {
  if((del3[i,4]-del3[i-1,4])<=6){
    del4 <- rbind(del4,ceiling(del3[i,]+del3[i-1,])/2);i=i+2
  }else {del4 <- rbind(del4,del3[i-1,]);i=i+1}
}
del4$bin <- ceiling(del4$bin)
del4[,5:29] <- ceiling(del4[,5:29])
fwrite(del4,"boundary_all_samples_filtered.bed",sep="\t",quote=F)
del4 <- cbind(del4,data.frame(rep(0,nrow(del4))))

#del5 <- del4[1,]
#del5 <- del5[-1,]
#i=2
#while (i <= nrow(del4)) {
#  if((del4[i,4]-del4[i-1,4])<=6){
#    del5 <- rbind(del5,ceiling(del4[i,]+del4[i-1,])/2);i=i+2
#  }else {del5 <- rbind(del5,del4[i-1,]);i=i+1}
#}
#del5$bin <- ceiling(del5$bin)
#del5[,5:29] <- ceiling(del5[,5:29])

stat <- t(data.frame(rep(0,25)))
rownames(stat) <- "num"
colnames(stat) <- c(1:25)
for (i in c(1:nrow(del4))){
  test <-  as.data.frame(table(del4[i,]>0))
  stat[1,test[2,2]-4]=stat[1,test[2,2]-4]+1
}
stat
png("TAD_Switch_all.png",width = 1100,height = 520)
barplot(c(stat/sum(stat)),col=cm.colors(30))
axis(1,seq(0.7,29.5,1.2),c(1:25),las=1,cex.lab=1)
dev.off()
