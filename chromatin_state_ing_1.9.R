##chromatin_state 多种来源基因组表观组类型对比
#input dir of files & dir of restored address


#install.packages("ggplot2")
library("ggplot2")
#install.packages('plotly')
library(plotly)
#install.packages("ggsci")
library("ggsci")

# Args <- commandArgs()
# myfiles <- Sys.glob(paste(Args[6],'/*',sep=''))
# len=length(myfiles)
# a=strsplit(myfiles,'[/]')
# filename=data.frame()
# for (i in c(1:len)){
#   c=a[[i]][length(a[[1]])]
#   filename[i,1]=c
# }

filltocolor <- c('#FF0000','#FF4500','#FF4500','#FF4500',
                 '#008000','#008000','#008000','#009600',
                 '#C2E105','#C2E105','#C2E105','#C2E105',
                 '#FFC34D','#FFC34D','#FFC34D','#FFFF00',
                 '#FFFF00','#FFFF00','#FFFF66','#66CDAA',
                 '#8A91D0','#E6B8B7','#7030A0','#808080',
                 '#FFFFFF')
for(i in c(1:nrow(filename))){
  #dat<-read.table(paste(Args[6],'/',filename[i,1],sep=''))
  dat<-read.table('/Users/guolu/Desktop/final_2_1e-6')
  df<-data.frame(
    SNP=dat$V3,  #????????????????????????????????????colname
    state=dat$V4, #state(cell??????????)
    segmentation=dat$V5,
    color_index=dat$V5
  )
  
  seg_1<-levels(factor(df$segmentation)) 
  color_index<-c()
  for (j in c(1:length(seg_1))) {
    color_index<-c(color_index,as.numeric(strsplit(seg_1[j],'[_]')[[1]][1]))
  }
  df$SNP<-as.character(df$SNP) #??????????????????
  
  p <- ggplot(df, aes(x=SNP, y=state)) +
    geom_tile(aes(fill = as.factor(df$segmentation)),color = 'white') + #color index have to conform with the level of segmentation
    scale_color_npg() +
    theme(
      #axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), #tick labels along axes
      axis.text.x = element_blank(),
      panel.background = element_blank(), #background underneath legend keys
      legend.key = element_blank() 
    ) + 
    #coord_equal(ratio = 1) + #y/x=1
    labs(x ="", y ="") +
    scale_fill_manual('STATE',values = filltocolor[color_index])+
    ggtitle('Chromatin State of Specific SNP')
  p<-p+theme(legend.key.size=unit(2,'mm'))
  p<-p+theme(legend.position='right',legend.key.width =unit(1,'mm'))
  
  #ggplotly(p) 报错
  #ggsave(paste(Args[7],'/',filename[i,1],'.svg',sep='')) 
  ggsave('/Users/guolu/Desktop/Chromatin State of Specific SNP.pdf')
}

# test_seg<-levels(factor(df$segmentation))
# b<-strsplit(test_seg,'[_]')
# d<-matrix(NA,25,2)
# for (i in 1:length(b)){
#   d[i,1]<-b[[i]][1]
#   d[i,2]<-b[[i]][2]
# }
# rownames(d)<-d[,1]
# d_1<-d[as.character(sort(as.numeric(d[,1]))),]
# d_2<-d_1[,2]
###
#Sys.setenv("PATH" = paste(Sys.getenv("PATH"),"/Applications/orca.app", sep = .Platform$path.sep))

