## analyse AB compartments difference across all samples after AB_Call.R run
## step1: plot all the samples
## step2: analyse

library(data.table)
setwd("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/matrix_40k/AB_Compartments")

### step1 plot all the samples
path <- "/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/matrix_40k/AB_Compartments/compartment_data/"
file <- list.files(path,pattern = "_compartments.txt")
filename <- paste0("CRC-",gsub("_CN_merged_40k_compartments.txt","",file[1]),"-N")
matrix <- fread(paste0(path,file[1]),select = c(3,4,5,8),col.names = c("chr","start","end",filename))
for (i in c(2:length(file))){
  filename <- substr(gsub("_merged_40k_compartments.txt","",file[i]),4,5)
  if (filename=="CN"){
    filename <- paste0("CRC-",substr(file[i],1,2),"-N")
  }else {if (filename=="CT"){
    filename <- paste0("CRC-",substr(file[i],1,2),"-T")
  }else {filename <- paste0("CRC-",substr(file[i],1,2),"-M")}}
  mid <- fread(paste0(path,file[i]),select = 8,col.names = filename)
  matrix <- cbind(matrix,mid)
  rm(mid)
  rm(filename)
}
mat_chr2 <- as.data.frame(matrix[matrix$chr=="chr2",])
mat_chr2[is.na(mat_chr2)] <- 0
mat_chr2 <- mat_chr2[1:2250,c(1:24)]
mat_chr2 <- mat_chr2[,c(1:4,6,9,12,15,20,22,5,7,10,13,16,18,19,21,23,8,11,14,17,24)]

#pdf("AB_chr2:1-90000000.pdf",width=1000,height = 800)
png("AB_chr2:1-90000000.png",width=1000,height = 800)
par(mfrow=c(21,1),mar=c(0,10,0.2,0),oma=c(0,0,0,0),mgp=c(1,0,0),xpd=NA)
for (i in c(4:ncol(mat_chr2))) {
  #max <- round(max(abs(mat_chr2[,i])),1)
  barplot(mat_chr2[,i],col = ifelse(mat_chr2[,i]>0, "darkred", "darkblue"),border = NA,axes = F,ylim = c(-0.02,0.02),
          space=0,cex.axis = 2)
  mtext(names(mat_chr2)[i],side = 2,cex = 1.5,las=2,line = -3)
#  axis(2,c(-0.02,0.02),c(-0.02,0.02),line = -2,las=2,cex.lab=1)
}
dev.off()

### step2 differential analysis
matrix <- matrix[,-c(1,2,3)]
matrix[is.na(matrix),] <- 0
del <- c()
k=1
for (i in c(1:nrow(matrix))){
  test <- as.data.frame(table(matrix[i,] == 0))
  condi <- nrow(test[test$Var1=="TRUE",])
  if(condi == 1){
    del[k] <- i
    k=k+1
  }
}
matrix <- matrix[-del,]

stat <- t(data.frame(rep(0,21)))
rownames(stat) <- "num"
#colnames(stat) <- c("25A","24A/1B","23A/2B","22A/3B","21A/4B","20A/5B","19A/6B","18A/7B","17A/8B","16A/9B","15A/10B","14A/11B",
#                    "13A/12B","12A/13B","11A/14B","10A/15B","9A/16B","8A/17B","7A/18B","6A/19B","5A/20B","4A/21B","3A/22B",
#                    "2A/23B","1A/24B","25B")
colnames(stat) <- c("20A","19A/1B","18A/2B","17A/3B","16A/4B","15A/5B","14A/6B","13A/7B","12A/8B","11A/9B","10A/10B",
                    "9A/11B","8A/12B","7A/13B","6A/14B","5A/15B","4A/16B","3A/17B","2A/18B","1A/19B","20B")
for (i in c(1:nrow(matrix))){
  test <-  as.data.frame(table(matrix[i,]>0))
  if(test[1,2]==20 & test[1,1]=="FALSE"){
    stat[1,1+test[1,2]]=stat[1,1+test[1,2]]+1
  }else stat[1,21-test[1,2]]=stat[1,21-test[1,2]]+1
}
pdf("AB_Switch_all.pdf",width = 10,height = 5)
color <- colorRampPalette(c("darkred","darkblue"))(21)
barplot(c(stat/sum(stat)),col=color,cex.axis = 1.5)
axis(1,seq(0.7,24.7,1.2),c("20A","19A/1B","18A/2B","17A/3B","16A/4B","15A/5B","14A/6B","13A/7B","12A/8B","11A/9B","10A/10B",
                           "9A/11B","8A/12B","7A/13B","6A/14B","5A/15B","4A/16B","3A/17B","2A/18B","1A/19B","20B")
                            ,las=1,cex.lab=1,cex.axis=1.5)
dev.off()
sum(stat[1,1],stat[1,21])/sum(stat)

### step3 analysis by tissues
func <- function(tissue){
  path <- "/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/matrix_40k/AB_Compartments/compartment_data/"
  file <- list.files(path,pattern = paste0(tissue,"_merged_40k_compartments.txt"))
  filename <- gsub("_merged_40k_compartments.txt","",file[1])
  matrix <- fread(paste0(path,file[1]),select = c(3,4,5,8),col.names = c("chr","start","end",filename))
  for (i in c(2:length(file))){
    filename <- gsub("_merged_40k_compartments.txt","",file[i])
    mid <- fread(paste0(path,file[i]),select = 8,col.names = filename)
    matrix <- cbind(matrix,mid)
    rm(mid)
    rm(filename)
  }
  n.sample <- length(file)
  matrix <- matrix[,-c(1,2,3)]
  matrix[is.na(matrix),] <- 0
  del <- c()
  k=1
  for (i in c(1:nrow(matrix))){
    test <- as.data.frame(table(matrix[i,] == 0))
    condi <- nrow(test[test$Var1=="TRUE",])
    if(condi == 1){
      del[k] <- i
      k=k+1
    }
  }
  matrix <- matrix[-del,]
  print(n.sample)
  return(matrix)
}

# plot CN
tissue <- "CN"
matrix <- lapply(tissue,func)
matrix <- as.data.frame(matrix)
stat <- t(data.frame(rep(0,7)))
rownames(stat) <- "num"
colnames(stat) <- c("6A","5A/1B","4A/2B","3A/3B","2A/4B","1A/5B","6B")
for (i in c(1:nrow(matrix))){
  test <-  as.data.frame(table(matrix[i,]>0))
  if(test[1,2]==6 & test[1,1]=="FALSE"){
    stat[1,1+test[1,2]]=stat[1,1+test[1,2]]+1
  }else stat[1,7-test[1,2]]=stat[1,7-test[1,2]]+1
}
pdf(paste0("AB_Switch_",tissue,".pdf"),width = 5,height = 5)
barplot(c(stat/sum(stat)),col=colorRampPalette(c("darkred","darkblue"))(9),cex.axis = 1.2)
axis(1,seq(0.7,7.9,1.2),c("6A","5A/1B","4A/2B","3A/3B","2A/4B","1A/5B","6B"),las=1,cex.axis=1.2)
#axis(1,seq(0.7,13.9,1.2),c("11A","10A/1B","9A/2B","8A/3B","7A/4B","6A/5B","5A/6B","4A/7B","3A/8B","2A/9B","1A/10B","11B"),las=1,cex.axis=1.2)
dev.off()
sum(stat[1,1],stat[1,7])/sum(stat)

# plot CT
tissue <- "CT"
matrix <- lapply(tissue,func)
matrix <- as.data.frame(matrix)
stat <- t(data.frame(rep(0,10)))
rownames(stat) <- "num"
colnames(stat) <- c("9A","8A/1B","7A/2B","6A/3B","5A/4B","4A/5B","3A/6B","2A/7B","1A/8B","9B")
for (i in c(1:nrow(matrix))){
  test <-  as.data.frame(table(matrix[i,]>0))
  if(test[1,2]==9 & test[1,1]=="FALSE"){
    stat[1,1+test[1,2]]=stat[1,1+test[1,2]]+1
  }else stat[1,10-test[1,2]]=stat[1,10-test[1,2]]+1
}
pdf(paste0("AB_Switch_",tissue,".pdf"),width = 7.5,height = 5)
barplot(c(stat/sum(stat)),col=colorRampPalette(c("darkred","darkblue"))(11),cex.axis = 1.2)
axis(1,seq(0.7,11.5,1.2),c("9A","8A/1B","7A/2B","6A/3B","5A/4B","4A/5B","3A/6B","2A/7B","1A/8B","9B"),las=1,cex.axis=1.2)
dev.off()
sum(stat[1,1],stat[1,10])/sum(stat)

# plot LT
tissue <- "LT"
matrix <- lapply(tissue,func)
matrix <- as.data.frame(matrix)
stat <- t(data.frame(rep(0,6)))
rownames(stat) <- "num"
colnames(stat) <- c("5A","4A/1B","3A/2B","2A/3B","1A/4B","5B")
for (i in c(1:nrow(matrix))){
  test <-  as.data.frame(table(matrix[i,]>0))
  if(test[1,2]==5 & test[1,1]=="FALSE"){
    stat[1,1+test[1,2]]=stat[1,1+test[1,2]]+1
  }else stat[1,6-test[1,2]]=stat[1,6-test[1,2]]+1
}
pdf(paste0("AB_Switch_",tissue,".pdf"),width = 5,height = 5)
barplot(c(stat/sum(stat)),col=colorRampPalette(c("darkred","darkblue"))(6),cex.axis = 1.2)
axis(1,seq(0.7,6.7,1.2),c("5A","4A/1B","3A/2B","2A/3B","1A/4B","5B"),las=1,cex.axis=1.2)
dev.off()
sum(stat[1,1],stat[1,6])/sum(stat)
