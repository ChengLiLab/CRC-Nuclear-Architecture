### plot heatmap of T vs N or M vs T

setwd("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/onTAD/N4T7_ontad")
options(scipen=200)
library(data.table)
library(pheatmap)

###T vs N
#ITGB1BP1 bou: chr2:9000001-10500000
#HOXD family: chr2:175400001-177400000
#RRM2 tad: chr2:10200000-10920000  
pat <- "05"
coef <- 1298.26/1318.22 #pat==01
coef <- 1284/1247.1 #pat==02
coef <- 1329.03/1342.29 #pat==03
coef <- 1219/1192.72 #pat==04
coef <- 1342.78/1292.74 #pat==05
coef <- 1261.84/1309.05 #pat==10

genename <- "MYC"
chr <- "chr8"
start <- 128200000
end <- 130500000
resolution <- 40000
col_cutoff <- 10
sample1 <- paste0(pat,"_CN")
matrix_CN <- as.data.frame(fread(paste0("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/down_sampling/TAD_50M/dense_matrix/headered/",sample1,"_40000_iced_",chr,"_normalmatrix.txt"),header=T))
matrix_CN <- matrix_CN[,-1]
matrix_CN <- matrix_CN*coef
#apply(matrix_CN[,100:105], 2,sum)
sample2 <- paste0(pat,"_CT")
matrix_CT <- as.data.frame(fread(paste0("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/down_sampling/TAD_50M/dense_matrix/headered/",sample2,"_40000_iced_",chr,"_normalmatrix.txt"),header=T))
matrix_CT <- matrix_CT[,-1]
#apply(matrix_CT[,100:105], 2,sum)

start = round(start/resolution)
end = round(end/resolution)
### extract Hi-C matrix and converted to triangle format
matrix_CN <- as.matrix(matrix_CN[(start+1):end,(start+1):end])
matrix_CN[matrix_CN>=col_cutoff]<-col_cutoff
matrix_CT <- as.matrix(matrix_CT[(start+1):end,(start+1):end])
matrix_CT[matrix_CT>=col_cutoff]<-col_cutoff
matrix <- matrix_CT-matrix_CN
range(matrix)
max <- ceiling(max(abs(min(matrix)),max(matrix)))
#breaks
bk <- c(seq(-max,-0.1,by=0.01),seq(0,max,by=0.01))
# plot heatmap
pdf(paste0(genename,"/",pat,"_compared_suqare.pdf"),3.2,3)
pheatmap(matrix,border_color = NA,cluster_rows = F,cluster_cols = F,show_colnames = F,show_rownames = F,
         color = c(colorRampPalette(colors = c("darkblue","white"))(length(bk)/2),colorRampPalette(colors = c("white","darkred"))(length(bk)/2)),
         legend_breaks=seq(-8,8,2),breaks=bk,main = paste0("CRC-",pat," (T-N)"),main.cex=1.3)
dev.off()

###M vs T
#FOXF1 gene TAD: chr16: 86000000-88000000
#RRAGD: chr6:90120000-90520000
pat <- "05"
coef <- 1263.84/1284 #pat==02
coef <- 1336.87/1329.03 #pat==03
coef <- 1272.05/1219 #pat==04
coef <- 1307.84/1342.78 #pat==05
coef <-  1332.72/1269.73 #pat==14

genename <- "RRAGD"
chr <- "chr6"
start <- 89000000
#start <- (round(start/40000)-40)*40000+1
end <- 91000000
#end <- (round(end/40000)+40)*40000
resolution <- 40000
col_cutoff <- 10
matrix <- as.data.frame(fread(paste0("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/down_sampling/TAD_50M/dense_matrix/headered/",pat,"_LT_40000_iced_",chr,"_normalmatrix.txt"),header=T))
matrix <- matrix[,-1]
matrix[matrix>=col_cutoff] <- col_cutoff
matrix2 <- as.data.frame(fread(paste0("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/down_sampling/TAD_50M/dense_matrix/headered/",pat,"_CT_40000_iced_",chr,"_normalmatrix.txt"),header=T))
matrix2 <- matrix2[,-1]
matrix2 <- coef*matrix2
matrix2[matrix2>=col_cutoff] <- col_cutoff
start = round(start/resolution)
end = round(end/resolution)
matrix <- as.matrix(matrix[(start+1):end,(start+1):end])
matrix2 <- as.matrix(matrix2[(start+1):end,(start+1):end])
matrix<- matrix-matrix2

range(matrix)
max <- ceiling(max(abs(min(matrix)),max(matrix)))
#breaks
bk <- c(seq(-max,-0.1,by=0.01),seq(0,max,by=0.01))
# plot heatmap
pdf(paste0(genename,"/",pat,"_compared_suqare.pdf"),3.2,3)
pheatmap(matrix,border_color = NA,cluster_rows = F,cluster_cols = F,show_colnames = F,show_rownames = F,
         color = c(colorRampPalette(colors = c("darkblue","white"))(length(bk)/2),colorRampPalette(colors = c("white","darkred"))(length(bk)/2)),
         legend_breaks=seq(-8,8,2),breaks=bk,main = paste0("CRC-",pat," (M-T)"),main.cex=1.3)
dev.off()

