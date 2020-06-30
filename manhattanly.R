#GWAS后续_manhattanly图
#input gwas.assoc.txt | anno.csv(注释文件)
#output manhattanly.html


install.packages("manhattanly")
if (!require("pacman")) install.packages("pacman")
pacman::p_load_gh("sahirbhatnagar/manhattanly")
library("manhattanly")
install.packages("tidyr")
library("tidyr")
install.packages("htmlwidgets")
library("htmlwidgets")

setwd("F:/2016_2018") 
myfiles <- Sys.glob("F:/2016_2018/gwas_result/*.assoc.txt")
len=length(myfiles)
a=strsplit(myfiles,'[.]')
filename=data.frame()
for (i in c(1:len)){
  b=strsplit(a[[i]][2],'[/]')
  n=length(b[[1]])
  c=b[[1]][n]
  filename[i,1]=c
}

anno<-read.table("F:/2016_2018/Uygur3M.anno.hg19_multianno.csv",sep=",",head=TRUE)
anno_1<-unite(anno, "Chr_Start", Chr,Start,remove = FALSE)

for (i in c(1:nrow(filename))){
  assoc_results <- read.table(paste("F:/2016_2018/gwas_result/Uygur3M.",filename[i,1],".assoc.txt",sep=''), head=TRUE)
  #È¡anno_assoc½»¼¯
  assoc_results_1<-unite(assoc_results, "chr_ps", chr,ps,remove = FALSE)
  anno_new<-anno_1[anno_1$Chr_Start %in% assoc_results_1$chr_ps,]
  #ºÏ²¢×¢ÊÍ£¬½«×¢ÊÍ¼ÓÈëassoc
  anno_new<-unite(anno_new,"annotation",Func.refGene,oreganno,Gene.refGene,wgEncodeRegTfbsClusteredV3)
  assoc_results_1$annotation<-anno_new$annotation
  #ÉèÖÃpãÐÖµ
  assoc_results_p5<-assoc_results_1[assoc_results_1$p_wald<1e-4,]
  #manhattanly ÐèÒªspecific name
  assoc_results_p5$CHR<-assoc_results_p5$chr
  assoc_results_p5$BP<-assoc_results_p5$ps
  assoc_results_p5$P<-assoc_results_p5$p_wald
  widget<-manhattanly(assoc_results_p5,snp="rs",gene="annotation")
  
  saveWidget(widget,file=paste("F:/2016_2018/feature manhattan/anno/",filename[i,1],".html",sep=""),selfcontained=TRUE,libdir="lib",background = "white", title=filename[i,1])
}
  
  
  
  