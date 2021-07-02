### normalized raw matrix at defined resolution, add header info and prepare input matrix for insulation socre

######define your parameters
options(scipen=200)
resolution <- 100000
human <- "data_01+02"
ID <- "HIC86T"

library(diffHic)
library(reshape2)
library(Homo.sapiens)
#library(BSgenome.Hsapiens.UCSC.hg19)
#seg.frags <- segmentGenome(BSgenome.Hsapiens.UCSC.hg19, size=40000) #segmentGenome function was abondand since 2017.10
human.genes = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
seg.frags = tileGenome(seqinfo(human.genes), tilewidth=resolution, cut.last.tile.in.chrom=T)

##########################
for(num in c(1:23)){
  if(num == 23){
    num = "X"
  }
  finder =  seg.frags[seqnames(seg.frags) == paste0("chr",num,collapse = "")]
  
  name = lapply("bin", paste,end(finder)/resolution,"|hg19|",as.character(seqnames(finder))[1],":",
                start(finder),"-",end(finder)+1,sep="")[[1]]
  name1 = name
  name1[1] = paste0("\t",name1[1],collapse = "")
  
  Matrix = fread(paste0("/lustre/user/liclab/ganjb/Projects/colon_cancer/",human,"/HiC/hicpro_results/hic_results/matrix/",ID,
                         "/iced/",resolution,"/dense_matrix/",ID,"_",resolution,"_iced_chr",num,"_dense.matrix"))
  #Matrix = Matrix[,-c(1:3)]
  
  rownames(Matrix) = name
  
  colnames(Matrix) = name1
  
  write.table(Matrix,file = paste0("/lustre/user/liclab/ganjb/Projects/colon_cancer/",human,"/HiC/hicpro_results/hic_results/matrix/",
                                   ID,"/iced/",resolution,"/dense_matrix/headered/",ID,"_",resolution,"_iced_chr",num,"_normalmatrix.txt",
                                   collapse = ""),quote = F,sep = "\t")
  
  print(num)
}
