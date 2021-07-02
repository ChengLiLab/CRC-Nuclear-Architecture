### draw Hi-C triangle heatmap and CNV and TAD boundaries

setwd("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/onTAD/N4T7_ontad")
options(scipen=200)
library(data.table)
#library(diffHic)
#library(BSgenome.Hsapiens.UCSC.hg19)

plotTAD <- function(sample,matrix,chr,start,end,resolution,boundary,ratio,col_cutoff) {
  start = round(start/resolution)
  end = round(end/resolution)
  
  ### extract Hi-C matrix and converted to triangle format
  raw_matrix <- as.matrix(matrix[(start+1):end,(start+1):end])
  raw_matrix[raw_matrix>=col_cutoff]<-col_cutoff
  vv = as.vector(raw_matrix)
  vv[is.na(vv)] = 0
  raw_matrix = matrix(vv,dim(raw_matrix))
  
  ti <- substr(sample,4,5)
  patient <- substr(sample,1,2)
  if (ti=="CN"){name <- paste0(patient,"-N")
  }else {if (ti=="CT"){name <- paste0(patient,"-T")}
    else {name <- paste0(patient,"-M")}}
  case_name <- paste0("CRC-",name,": ",chr,":",start*resolution/1000000,"mb-",end*resolution/1000000,"mb")
  
  re_matrix = array(NA,dim = dim(raw_matrix))
  mm = ll = dim(re_matrix)[1]
  k=0
  nn = 0
  for(j in 1:ll){
    if(j %% 2 == 0){k = k+1}
    for(i in 1:mm){
      re_matrix[ll-j+1,i+k] = raw_matrix[i,i+nn]
    }
    mm = mm -1
    nn = nn +1
  }
  
  ##################
  # start plotting 
  ##################
  nf <- layout(matrix(c(1,2,3),3,1,byrow = TRUE),heights = c(4,0.3,1.8))
  
  #  plot triangle heatmap
  par(mar=c(0,6,3,2))
  row_num <- dim(re_matrix)[1]
  col_num <- dim(re_matrix)[2]
  
  re = as.vector(t(re_matrix))
  re = re[is.na(re) == F]
  #generate color vector
  color <- colorRampPalette(c("white","red"))
  re <- re*1000+1
  color <- color(col_cutoff*1000+1)[round(re)]
  source("/lustre/user/liclab/ganjb/scripts/plot_triangle_heatmap_function.R")
  m =dim(re_matrix)[1]
  exam_52 <- construct_xy(m)
  
  plot(x=-col_num:col_num,y=seq(0,m,by=0.5),type="n",frame.plot = F,cex.axis=1,axes = FALSE,
       cex.lab=1.6,ylab = "Contact matrix",xlab = "",main = case_name,cex.main = 2)
  #nn = round((end - start)/5)
  #axis(1,at = seq(-col_num,col_num,by = 2*nn),labels = unlist(lapply(prettyNum((seq(start,end,by = nn))*resolution/1000000,big.mark=","),paste0,"mb")),
  #     cex.axis=1.3)
  polygon(exam_52[,1],exam_52[,2],col = color,border = NA) 
  
  ### plot TAD boundaries
  par(mar=c(0,6,0,2))
 # plot_TAD_boundary <- function(boundary,sample,chr,start,end){
    library(data.table)
    library(plyr)
    options(scipen = 200)
    #start plotting
    aa <- (boundary$chr==chr)
    chr7 <- as.data.frame(boundary[aa,])
    bb <- (chr7$start >= (start*resolution-1))
    chr7 <- chr7[bb,]
    cc <- (chr7$end <= end*resolution)
    chr7 <- chr7[cc,]
    chr7[chr7 >= 1] <- 1
    i <- (colnames(chr7)==sample)
    plot(x=1:nrow(chr7),y=seq(length=nrow(chr7),from=0,to=1),type="n",frame.plot = F,cex.axis=1,axes = FALSE,
         cex.lab=1.6,ylab = "",xlab = "",main = "",cex.main = 2)
    barplot(chr7[,i],col = "darkblue",border = ifelse(chr7[,i]==1,"darkblue","white"),axes = F,ylim = c(0.6,1),space=0,cex.axis = 2,
            cex.lab=1.5,xlab= "",ylab = "TAD",add = T)
    #mtext(names(chr7)[i],side = 2,cex = 1,las=2,line=-3)
 # }
    
  ### plot CNV
  par(mar=c(4,6,0,2))
  i <- strsplit(chr,"chr")[[1]][2]
  ploidy = 2
  maxLevelToPlot = 2
  tt <- which(ratio$Chromosome==i)
  tt = tt[(start+1):(end)]
  cnv <- ratio[tt,]
  if (length(tt)>0) {
    cnv$Start <- (cnv$Start-1)/resolution
    plot(cnv$Start,cnv$Ratio*ploidy,ylim = c(0,maxLevelToPlot*ploidy),frame.plot = F,cex.lab=1.5,
         xlab = "",axes = F,ylab = "CNV",cex.lab=1.5,pch = ".",col = colors()[88],cex = 5)
    nn = round((end - start)/5)
    axis(1,at = seq(cnv$Start[1],(cnv$Start[nrow(cnv)]+1),length.out = 5),labels = unlist(lapply(prettyNum((seq(start,end,length.out = 5))*resolution/1000000,big.mark=","),paste0,"mb")))
    #axis(1,labels = T)
    axis(2,at=c(0,2,4),labels = c(0,2,4))
    tt <- which(cnv$Chromosome==i  & cnv$CopyNumber>ploidy )
    points(cnv$Start[tt],cnv$Ratio[tt]*ploidy,pch = ".",col = colors()[136],cex = 5)
    #text(seq(0,col_num*(10^8)/100,by=(10^7)),-1.5,prettyNum(seq(0,col_num*(10^2)/100,by=10),big.mark=","),srt=-45,xpd = TRUE)
    
    tt <- which(cnv$Chromosome==i  & cnv$cnv==maxLevelToPlot & cnv$CopyNumber>ploidy)	
    points(cnv$Start[tt],cnv$Ratio[tt]*ploidy,pch = ".",col = colors()[136],cex=5)
    
    tt <- which(cnv$Chromosome==i  & cnv$CopyNumber<ploidy & cnv$CopyNumber!= -1)
    points(cnv$Start[tt],cnv$Ratio[tt]*ploidy,pch = ".",col = colors()[461],cex = 5)
    tt <- which(cnv$Chromosome==i)
    #UNCOMMENT HERE TO SEE THE PREDICTED COPY NUMBER LEVEL:
    points(cnv$Start[tt],cnv$CopyNumber[tt], pch = ".", col = colors()[240],cex=3)
    #nn = round((end - start)/5)
    #axis(1,at = seq(1,nrow(cnv),by = nn),labels = unlist(lapply(prettyNum((seq(start,end,by = nn))*resolution/1000000,big.mark=","),paste0,"mb")),
    #     cex.axis=1.1)
  }
}

#ITGB1BP1 bou: chr2:7520001-11560000
#HOXD family TAD: chr2:176105001-177030000
#HOXD family TAD: chr2:176160000-177000000
#FOXF1 gene TAD: chr16: 86200000-86600000
#CCL5 gene: chr17:34198495-34207377
#CCL22 gene: chr16:57392694-57400102
#IFIT2 gene: chr10:91061705-91069033
#BMP5: chr6:55618450-55740388
#RRAGD tad: chr6:90120000-90520000

sample <- "04_CT"
chr <- "chr6"
start <- 90120000       
start <- (round(start/40000)-10)*40000+1
end <- 90520000
end <- (round(end/40000)+10)*40000
resolution <- 40000
col_cutoff <- 10
matrix <- as.data.frame(fread(paste0("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/down_sampling/TAD_50M/dense_matrix/headered/",sample,"_40000_iced_",chr,"_normalmatrix.txt"),header=T))
matrix <- matrix[,-1]
ratio <-read.table(paste0("/lustre/user/liclab/ganjb/Projects/colon_cancer/CNV/CNV_data/",sample,"_ratio.txt"),header=TRUE)
#boundary <- fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/onTAD/colon_ontad_2.1/12tad_boundary_level_3bin_50/boundary_all_samples.txt")
boundary <- fread("/lustre/user/liclab/ganjb/Projects/colon_cancer/HiC/down_sampling/TAD_50M/TAD_analysis/boundary_all_samples.bed")
plotTAD(sample,matrix,chr,start,end,resolution,boundary,ratio,col_cutoff)

