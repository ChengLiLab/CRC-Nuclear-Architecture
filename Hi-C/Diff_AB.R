## analyse AB compartments different after AB_Call.R run
## step1: differential analysis between samples
## step2: plot

library(data.table)
setwd("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/merged_matrix_40k/AB_Compartments")

### step1 find difference between samples
CN <- fread("../CN_merged_40k_compartments.txt",select = c(3,4,5,10),col.names = c("chr","start","end","CN"))
CT <- fread("../CT_merged_40k_compartments.txt",select = 10,col.names = c("CT"))
LT <- fread("../LT_merged_40k_compartments.txt",select = 10,col.names = c("LT"))

mat <- cbind(cbind(CN,CT),LT)
A2B_NT <- mat[mat$CN=="A" & mat$CT=="B"]
B2A_NT <- mat[mat$CN=="B" & mat$CT=="A"]
A2B_TL <- mat[mat$CT=="A" & mat$LT=="B"]
B2A_TL <- mat[mat$CT=="B" & mat$LT=="A"]
A2B2A <- A2B_NT[A2B_NT$LT=="A"]
B2A2B <- B2A_NT[B2A_NT$LT=="B"]
fwrite(A2B_NT,"A2B_NT",quote = F)
fwrite(B2A_NT,"B2A_NT",quote = F)
fwrite(A2B_TL,"A2B_TL",quote = F)
fwrite(B2A_TL,"B2A_TL",quote = F)
fwrite(A2B2A,"A2B2A",quote = F)
fwrite(B2A2B,"B2A2B",quote = F)
#nrow(A2B2A)==nrow(B2A_TL[B2A_TL$CN=="A"])
#nrow(B2A2B)==nrow(A2B_TL[A2B_TL$CN=="B"])
AA_NT_len <- nrow(mat[mat$CN=="A" & mat$CT=="A"])
BB_NT_len <- nrow(mat[mat$CN=="B" & mat$CT=="B"])
AA_TL_len <- nrow(mat[mat$CT=="A" & mat$LT=="A"])
BB_TL_len <- nrow(mat[mat$CT=="B" & mat$LT=="B"])

A2B_NT_len <- nrow(A2B_NT)
B2A_NT_len <- nrow(B2A_NT)
A2B_TL_len <- nrow(A2B_TL)
B2A_TL_len <- nrow(B2A_TL)

### plot AB compartments of 3 samples
AAA_len <- nrow(mat[mat$CN=="A" & mat$CT=="A" & mat$LT=="A"])
BBB_len <- nrow(mat[mat$CN=="B" & mat$CT=="B" & mat$LT=="B"])
AAB_len <- nrow(mat[mat$CN=="A" & mat$CT=="A" & mat$LT=="B"])
ABA_len <- nrow(mat[mat$CN=="A" & mat$CT=="B" & mat$LT=="A"])
ABB_len <- nrow(mat[mat$CN=="A" & mat$CT=="B" & mat$LT=="B"])
BBA_len <- nrow(mat[mat$CN=="B" & mat$CT=="B" & mat$LT=="A"])
BAA_len <- nrow(mat[mat$CN=="B" & mat$CT=="A" & mat$LT=="A"])
BAB_len <- nrow(mat[mat$CN=="B" & mat$CT=="A" & mat$LT=="B"])
stat <- data.frame(CN=c(BBB_len,300,BBA_len,300,BAB_len,300,BAA_len,300,ABB_len,300,ABA_len,300,AAB_len,300,AAA_len),
                   CT=c(BBB_len,300,BBA_len,300,BAB_len,300,BAA_len,300,ABB_len,300,ABA_len,300,AAB_len,300,AAA_len),
                   LT=c(BBB_len,300,BBA_len,300,BAB_len,300,BAA_len,300,ABB_len,300,ABA_len,300,AAB_len,300,AAA_len))
stat <- as.matrix(stat)
par(mfrow=c(1,3),oma=c(3,5,0,5),mgp=c(3,2,0),xpd=T)
barplot(matrix(stat[,1]),col = c("darkblue","white","darkblue","white","darkblue","white","darkblue","white","darkred","white","darkred","white","darkred","white","darkred"),border=NA,axes = F)
text(0.75,0,"CN",pos=1,cex = 2.5)
legend(0,40000,c("A","B"),col=c("darkred","darkblue"),bty="n", pch=15,cex = 2.8,xjust = .55,xpd=T)
barplot(matrix(stat[,2]),col = c("darkblue","white","darkblue","white","darkred","white","darkred","white","darkblue","white","darkblue","white","darkred","white","darkred"),border=NA,axes = F)
text(0.75,0,"CT",pos=1,cex = 2.5,xpd=T)
#text(0.75,68000,"AB Comopartments Switch",pos = 3,cex = 4)
barplot(matrix(stat[,2]),col = c("darkblue","white","darkred","white","darkblue","white","darkred","white","darkblue","white","darkred","white","darkblue","white","darkred"),border=NA,axes = F)
text(0.75,0,"LT",pos=1,cex = 2.5)
text(1.4,36000,"7.912%",pos=1,cex = 1.5)

### plot AB compartments of 2 paired samples
stat <- data.frame(CN=c(BB_NT_len,300,B2A_NT_len,300,A2B_NT_len,300,AA_NT_len),
                   CT=c(BB_NT_len,300,B2A_NT_len,300,A2B_NT_len,300,AA_NT_len))
par(mfrow=c(1,2),oma=c(3,5,0,5),mgp=c(3,2,0),xpd=T)
barplot(matrix(stat[,1]),col = c("darkblue","white","darkblue","white","darkred","white","darkred"),axes=F,border=NA)
text(0.75,0,"CN",pos=1,cex = 2.5)
barplot(matrix(stat[,1]),col = c("darkblue","white","darkred","white","darkblue","white","darkred"),axes=F,border=NA)
text(0.75,0,"CT",pos=1,cex = 2.5,xpd=T)
#A2B_NT_len/sum(AA_NT_len,BB_NT_len,A2B_NT_len,B2A_NT_len)
#B2A_NT_len/sum(AA_NT_len,BB_NT_len,A2B_NT_len,B2A_NT_len)


stat <- data.frame(CT=c(BB_TL_len,300,B2A_TL_len,300,A2B_TL_len,300,AA_TL_len),
                   LT=c(BB_TL_len,300,B2A_TL_len,300,A2B_TL_len,300,AA_TL_len))
par(mfrow=c(1,2),oma=c(3,5,0,5),mgp=c(3,2,0),xpd=T)
barplot(matrix(stat[,1]),col = c("darkblue","white","darkblue","white","darkred","white","darkred"),axes=F,border=NA)
text(0.75,0,"CT",pos=1,cex = 2.5)
barplot(matrix(stat[,1]),col = c("darkblue","white","darkred","white","darkblue","white","darkred"),axes=F,border=NA)
text(0.75,0,"LT",pos=1,cex = 2.5,xpd=T)
#A2B_TL_len/sum(AA_TL_len,BB_TL_len,A2B_TL_len,B2A_TL_len)
#B2A_TL_len/sum(AA_TL_len,BB_TL_len,A2B_TL_len,B2A_TL_len)
