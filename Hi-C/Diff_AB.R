### analyse the difference of AB compartments between N and T or T and M 
### plot pie plot and A/B heatmap

## step1: merge all samples 
## step2: differential analysis
## step3: plotting

library(data.table)
setwd("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/matrix_40k/AB_Compartments")

###merge all the samples
path <- "/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/matrix_40k/AB_Compartments/"
file <- list.files(path,pattern = "_compartments.txt")
filename <- paste0("CRC",gsub("_CN_merged_40k_compartments.txt","",file[1]),"-N")
matrix <- fread(file[1],select = c(3,4,5,8),col.names = c("chr","start","end",filename))
for (i in c(2:length(file))){
  filename <- substr(gsub("_merged_40k_compartments.txt","",file[i]),4,5)
  if (filename=="CN"){
    filename <- paste0("CRC",substr(file[i],1,2),"-N")
  }else {if (filename=="CT"){
    filename <- paste0("CRC",substr(file[i],1,2),"-T")
  }else {filename <- paste0("CRC",substr(file[i],1,2),"-M")}}
  mid <- fread(file[i],select = 8,col.names = filename)
  matrix <- cbind(matrix,mid)
  rm(mid)
  rm(filename)
}
matrix <- as.data.frame(matrix)
rownames(matrix) <- paste0(matrix$chr,"_",matrix$start,"_",matrix$end)
matrix <- matrix[,-c(1:3)]
fwrite(matrix,"all_sample_40k_compartments_PC1.txt",col.names=T,row.names=T,quote=)
matrix[is.na(matrix)] <- 0
matrix[matrix<0] <- -1
matrix[matrix>0] <- 1
fwrite(matrix,"all_sample_40k_compartments.txt",col.names=T,row.names=T,quote=F)
matrix$N <- apply(matrix[,c(1,3,6,9,12,17,19)],1,sum)
matrix$T <- apply(matrix[,c(2,4,7,10,13,15,16,18,20)],1,sum)
matrix$M <- apply(matrix[,c(5,8,11,14,21)],1,sum)

A2B_NT <- matrix[(matrix$T-matrix$N)<=-8,]
A2B_NT <- separate(data.frame(name=rownames(A2B_NT)),"name",c("chr","start","end"),"_")
B2A_NT <- matrix[(matrix$T-matrix$N)>=8,]
B2A_NT <- separate(data.frame(name=rownames(B2A_NT)),"name",c("chr","start","end"),"_")

matrix$T <- apply(matrix[,c(4,7,10,13,20)],1,sum)
matrix$M <- apply(matrix[,c(5,8,11,14)],1,sum)
A2B_TM <- matrix[(matrix$M-matrix$T)<=-7,]
A2B_TM <- separate(data.frame(name=rownames(A2B_TM)),"name",c("chr","start","end"),"_")
B2A_TM <- matrix[(matrix$M-matrix$T)>=7,]
B2A_TM <- separate(data.frame(name=rownames(B2A_TM)),"name",c("chr","start","end"),"_")
fwrite(A2B_NT,"./diff/A2B_NT.bed",sep="\t",col.names =F,quote=F)
fwrite(B2A_NT,"./diff/B2A_NT.bed",sep="\t",col.names =F,quote=F)
fwrite(A2B_TM,"./20200909/A2B_TM.bed",sep="\t",col.names =F,quote=F)
fwrite(B2A_TM,"./20200909/B2A_TM.bed",sep="\t",col.names =F,quote=F)

###pie plot
NT <- c(75918-nrow(A2B_NT)-nrow(B2A_NT),nrow(A2B_NT),nrow(B2A_NT))
NT <- round(100*NT/75918,2)
NT_lab <- c("No Switch: ","A2B: ","B2A: ")
NT_lab <- paste0(NT_lab,NT,"%")
NT <- 100*NT
pie(NT,NT_lab,col=c("grey","darkblue","darkred"),border = NA)
TM <- c(75918-nrow(A2B_TM)-nrow(B2A_TM),nrow(A2B_TM),nrow(B2A_TM))
TM <- round(100*TM/75918,2)
TM_lab <- c("No Switch: ","A2B: ","B2A: ")
TM_lab <- paste0(TM_lab,TM,"%")
TM <- 100*TM
pie(TM,TM_lab,col=c("grey","darkblue","darkred"),border = NA)

###heatmap
library(pheatmap)
matrix <- as.data.frame(fread("all_sample_40k_compartments_PC1.txt"))
row.names(matrix) <- matrix$V1
matrix <- matrix[,-1]
matrix[is.na(matrix)==T] <- 0
region <- matrix[which(rownames(matrix)=="chr1_152080001_152120000"):which(rownames(matrix)=="chr1_165000001_165040000"),
                 c("CRC01-N","CRC02-N","CRC03-N","CRC04-N","CRC05-N","CRC10-N","CRC14-N","CRC01-T","CRC02-T","CRC03-T",
                   "CRC04-T","CRC05-T","CRC06-T","CRC07-T","CRC10-T","CRC14-T","CRC02-M","CRC03-M","CRC04-M","CRC05-M","CRC14-M")]
colnames(region) <- c("CRC-01-N","CRC-02-N","CRC-03-N","CRC-04-N","CRC-05-N","CRC-10-N","CRC-14-N","CRC-01-T","CR-C02-T",
                      "CRC-03-T","CRC-04-T","CRC-05-T","CRC-06-T","CRC-07-T","CRC-10-T","CRC-14-T","CRC-02-M","CRC-03-M",
                      "CRC-04-M","CRC-05-M","CRC-14-M")
#region <- t(scale(region,center = F,scale = T))
range(region)
region[region>=0.045] <- 0.045; region[region<=-0.045] <- -0.045
pheatmap(t(region),col = colorRampPalette(c("darkblue","white", "darkred"))(100),border_color = NA,show_colnames = F,
         cluster_rows = F,cluster_cols = F)
