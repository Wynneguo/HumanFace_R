#Manhattan图的两种画法
#input assoc.txt
#output gwas.pdf

##Manhattan_single sample

myfiles <- Sys.glob("F:/2016_2018/gwas_result/Uygur3M.6.assoc.txt")
#myfiles <- Sys.glob("/assoc_address/*.assoc.txt")
len=length(myfiles)
a=strsplit(myfiles,'[/]')
filename=data.frame()
for (i in c(1:len)){
  c=a[[1]][4]
  filename[i,1]=c
}
library("qqman")
color_set <- c("#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#66C2A5","#FC8D62")
for (i in c(1:nrow(filename))){
  results_log <- read.table('F:/2016_2018/gwas_result/Uygur3M.6.assoc.txt', head=TRUE)
  jpeg(paste("C:/Users/GuoLu/Desktop/",filename[i,1],"_pwald.jpeg",sep=''))
  #manhattan(results_log,chr="chr",bp="ps",p="p_wald",snp="rs",col = c("green4", "orange3"), main = paste("Manhattan plot: Uygur3M",filename[i,1],sep=''),ylim=c(0,10), chrlabs = as.character(c(1:22)), cex = 0.6, cex.axis = 0.9)
  manhattan(results_log,chr="chr",bp="ps",p="p_wald",snp="rs",col=color_set, main = paste("Manhattan plot: Uygur3M",filename[i,1],sep=''),ylim=c(0,10), cex = 0.6, cex.axis = 0.9,chrlabs = as.character(c(1:22)))
  
  dev.off()
}

##Manhattan_all_samples
myfiles <- Sys.glob("F:/2016_2018/gwas_result/*.txt")
#myfiles <- Sys.glob("/assoc_address/*.assoc.txt")
len=length(myfiles)
a=strsplit(myfiles,'[/]')
filename=data.frame()
for (i in c(1:len)){
  c=a[[i]][4]
  filename[i,1]=c
}

jpeg(paste("C:/Users/GuoLu/Desktop/all_sample_pwald.jpeg",sep=''))
for (i in c(1:nrow(filename))){
  results_log <- read.table(paste('F:/2016_2018/gwas_result/',filename[i,1],sep=''), head=TRUE)
  #manhattan(results_log,chr="chr",bp="ps",p="p_wald",snp="rs",col = c("green4", "orange3"), main = paste("Manhattan plot: Uygur3M",filename[i,1],sep=''),ylim=c(0,10), chrlabs = as.character(c(1:22)), cex = 0.6, cex.axis = 0.9)
  manhattan(results_log,chr="chr",bp="ps",p="p_wald",snp="rs",col=color_set, main = paste("Manhattan plot: ",filename[i,1],sep=''),ylim=c(0,10), cex = 0.6, cex.axis = 0.9,chrlabs = as.character(c(1:22)))
  par(new=T)
}
dev.off()


##ggplot for manhattan plot
#test pass
library('doBy')
gwasResults<-results_log
gwasResults$SNP1 <- seq(1, nrow(gwasResults), 1) 
gwasResults$chr <- factor(gwasResults$chr, levels = unique(gwasResults$chr))
chr <- summaryBy(SNP1~chr, gwasResults, FUN = median)
#ggplot2 ×÷Í¼
p <- ggplot(gwasResults, aes(SNP1, -log(p_wald, 10), color = chr)) +
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), panel.background = element_rect(fill = 'transparent')) +
  geom_point(aes(color = chr), show.legend = FALSE) +
  labs(x = 'Chromosome', y = expression(''~-log[10]~'(p_wald)')) +
  scale_x_continuous(breaks = chr$SNP1.median, labels = chr$chr, expand = c(0, 0)) +
  geom_hline(yintercept = c(5, 7), color = c('blue', 'red'), size = 0.5)
ggsave('C:/Users/GuoLu/Desktop/GWAS.pdf', p, width = 10, height = 4)
