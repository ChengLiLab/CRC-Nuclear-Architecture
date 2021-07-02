### plot heatmap for single sample

setwd("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/onTAD/N4T7_ontad")
options(scipen=200)
library(data.table)
library(pheatmap)

###T vs N
#RRM2 tad(up): chr2:10200000-10920000 (plot:chr2:9500000-11000000)
#ENTPD5 tad(down): chr14:74120000-74720000 (plot:chr14:73400000-74800000)
pat <- "05"
sample1 <- paste0(pat,"_CN")
genename <- "MYC"
chr <- "chr8"
start <- 128200000
end <- 130500000
resolution <- 40000
#col_cutoff <- 10
matrix_CN <- as.data.frame(fread(paste0("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/down_sampling/TAD_50M/dense_matrix/headered/",sample1,"_40000_iced_",chr,"_normalmatrix.txt"),header=T))
matrix_CN <- matrix_CN[,-1]

start = round(start/resolution)
end = round(end/resolution)

### extract Hi-C matrix
matrix_CN <- as.matrix(matrix_CN[(start+1):end,(start+1):end])

matrix <- matrix_CN
range(matrix)
col_cutoff <- as.numeric(quantile(matrix,0.9))
matrix[matrix>col_cutoff] <- col_cutoff
max <- ceiling(max(abs(min(matrix)),max(matrix)))
#breaks
bk <- c(seq(0,max,by=0.01))
# plot heatmap
pdf(paste0(genename,"/",sample1,"_single_suqare.pdf"),3.2,3)
pheatmap(matrix,border_color = NA,cluster_rows = F,cluster_cols = F,show_colnames = F,show_rownames = F,
         color = c(colorRampPalette(colors = c("white","darkred"))(length(bk))),
         legend_breaks=seq(0,max,3),breaks=bk,main = paste0(sample1),main.cex=1.3)
dev.off()
