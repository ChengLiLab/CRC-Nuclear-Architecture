### Used for plot TADs number and TADs lenght distribution

library(data.table)
library(parallel)
library(plyr)
options(scipen = 200)
setwd("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/matrix_40k/TAD_plot")

## process samples
func <- function(sample){
  library(data.table)
  all <- matrix(NA,0,1)
  out <- matrix(1:23,23,2,byrow = F)
  path <- "/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/matrix_40k/TAD/test/"
  file <- list.files(path,pattern = "boundaries.bed")
  file <- file[substring(file,4,8)==sample]
  file <- file[substring(file,21,24)!="chrY"]
  for(num in 1:length(file)){
    mat <- fread(paste0(path,"/",file[num]), header = F)
    colnames(mat) <- c("chrom","start","end","header","score")
    bou <- matrix((mat$start+mat$end)/2)
    len <- matrix(NA,nrow(bou)-1,1)
    len <- as.data.frame(cbind(len,mat[2:nrow(mat),4]))
    
    for(i in 2:nrow(bou)){    
      len[i-1,1] <- (bou[i,1]-bou[i-1,1])
    }
    len <- len[is.na(len$V1)==F,]
    #len <- len[len$V1<=3000000,] ##Delete TADs whose length > 2500000
    #len <- len[len$V1>=200000,]  ##Delete TADs whose length < 200000
    length(which(len<200000)) 
    out[num,2] <- nrow(len)
    all <- rbind(all,len)
  }
  print(sum(out[,2]))
  colnames(all) <- c("TAD_length","TAD_boundary")
  #fwrite(all,paste0(sample,"_filtered_boundaries.bed"),sep = "\t")
  #fwrite(out,paste0(sample,"_TAD_num.txt"),sep = "\t",col.names = F)
  return(all[,1])
}
cl <- makeCluster(9)
#samples <- c("01_CN","01_CT","02_CN","02_CT","02_LT","02_LN","03_CN","03_CT","03_LT","03_NT","04_CN","04_CT","04_LT","05_CN",
#             "05_CT","05_LT","06_CN","06_CT","07_CN","07_CT","08_CN","08_CT","09_CN","09_CT","10_CN","10_CT","11_CN","11_CT",
#             "13_CN","13_CT","13_LT","14_CN","14_CT","14_LT")
samples <- c("01_CN","01_CT","02_CN","02_CT","02_LT","03_CN","03_CT","03_LT","04_CN","04_CT","04_LT","05_CN",
             "05_CT","05_LT","06_CT","07_CT","10_CN","10_CT","11_CN","11_CT","14_CN","14_CT","14_LT")
res <- parLapply(cl,samples,func)
stopCluster(cl)

## plot TAD length density distribution and number
names(res) <- c("CN","CT","CN","CT","LT","CN","CT","LT","CN","CT","LT","CN",
                "CT","LT","CT","CT","CN","CT","CN","CT","CN","CT","LT")
col<-c()
CN <- data.frame(CN=c())
CT <- data.frame(CT=c())
LT <- data.frame(LT=c())
for (i in 1:length(names(res))) {
  if(names(res[i])=="CN"){
    col[i]="lightblue";CN <- rbind(CN,res[i])
    }else if(names(res[i])=="CT"){
      col[i]="orange";CT <- rbind(CT,res[i])
      }else if(names(res[i])=="LT"){
        col[i]="pink";LT <- rbind(LT,res[i])
        }else col[i]="grey"
}
colnames(CN) <- "N"
colnames(CT) <- "T"
colnames(LT) <- "M"
density(CN$N)
density(CT$T)
density(LT$M)
#png("TAD_length_distribution.png",width = 880,height = 510)
#par(mar=c(6, 6, 3.2, 3))
#plot(density(res[[1]]), main = "Distribution of TAD Length in All Samples",cex.main=2.5,xlab = "TAD Length",cex.lab=2,
#     col=col[1],xlim=c(0,3000000),ylim=c(0,0.0000012),cex.axis=1.5)
#for (i in 2:length(res)) {
#  lines(density(res[[i]]),col=col[i])
#}
#legend(2000000,0.0000012,c("N: 10.88M","T: 10.64M","M: 11.04M"),cex=2,col=c("lightblue","orange","pink"),lty=c(1,1),bty = "n")
#dev.off()

#### plot boxplot of TAD number
boxplot(list("N"=CN$N,"T"=CT$T,"M"=LT$M),col=c("lightblue","orange","pink"),
        ylab="length",cex.lab=1.5,cex.axis=1.2,outline = F,varwidth = T)

TAD_num <- t(data.frame(rep(0,23)))
colnames(TAD_num) <- c("N","T","N","T","M","N","T","M","N","T","M","N",
                       "T","M","T","T","N","T","N","T","N","T","M")
for (i in c(1:length(res))){
  TAD_num[1,i] <- length(res[[i]])
}
data <- data.frame(Number=TAD_num[,1:23],Tissue=colnames(TAD_num))
x1 = factor(data$Tissue, levels=c("N","T","M"))
boxplot( Number ~ x1, data, varwidth = T,col=c("lightblue","orange","pink"),
         xlab="Tissue", ylab="Number",cex.lab=1.5,cex.axis=1.2)

###plot conservation of TAD and AB
AB <- fread("./AB_conserve.txt")
AB$V2 <- 100-AB$V2
boxplot(list("Non-Met"=AB$V2[1:2],"Met"=AB$V2[3:7]),col=c("lightgrey","darkgrey"),ylab="Compartment Conservasion")
TAD <- fread("./TAD_conserve.txt")
boxplot(list("Non-Met"=TAD$T[1:2],"Met"=TAD$T[3:7]),col=c("lightgrey","darkgrey"),ylab="TAD Conservasion")

