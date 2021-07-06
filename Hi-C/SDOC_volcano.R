###plot volcano plot for SDOC

setwd("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/SDOC")
library(data.table)

all <- fread("NT_SDOC_P_value.csv")
all$p_value2 <- -log10(all$p_value)
up <- all[all$p_value<=0.05 & all$mean_T_N>0,]
down <- all[all$p_value<=0.05 & all$mean_T_N<0,]
plot(all$mean_T_N,all$p_value2,pch=20,col="black",lwd=1,ylim=c(0,2),
     xlab="mean(T - N)",ylab="-log10(p_Value)",cex.lab =1.5,asp=0,cex.axis=1.3,bty="l")
points(up$mean_T_N,up$p_value2,pch=20,col="orange",asp=0,cex=1)
points(down$mean_T_N,down$p_value2,pch=20,col="lightblue",asp=0,cex=1)

all <- fread("TM_SDOC_P_value.csv")
all$p_value2 <- -log10(all$p_value)
up <- all[all$p_value<=0.05 & all$mean_M_T>0,]
down <- all[all$p_value<=0.05 & all$mean_M_T<0,]

plot(all$mean_M_T,all$p_value2,pch=20,col="black",lwd=1,ylim=c(0,2),
     xlab="mean(M - T)",ylab="-log10(p_Value)",cex.lab =1.5,asp=0,cex.axis=1.3,bty="l")
points(up$mean_M_T,up$p_value2,pch=20,col="pink",asp=0,cex=1)
points(down$mean_M_T,down$p_value2,pch=20,col="orange",asp=0,cex=1)
