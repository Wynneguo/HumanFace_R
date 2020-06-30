#GWAS后分析
#input GWAS结果
#output  显著SNPsig_snps.txt ｜ LD_file | Locuszoom_file | Paintor File


library(data.table)
plink_path <- '/asnas/chenhua_group/Software/plink-1.9-x86_64/plink'
locuszoom_path <- '/asnas/chenhua_group/dongyi/SOFTWARES/locuszoom/bin/locuszoom' #´ý¸üÐÂÎªguoluÂ·¾¶
Paintor_path <- '/asnas/chenhua_group/guolu/toolkit/PAINTOR_V3.0-master/PAINTOR'
Uygur3M_path <- '/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Uygur3M'



myfiles <- Sys.glob("/asnas/chenhua_group/guolu/2016_2018/dongyi_data/GWAS_result/bolt_assoc_resultprocess1/*")
myfiles_signal<-character()
for (i in 1:length(myfiles)){
  stats <- fread(myfiles[i],header = T)
  stats <- stats[order(stats$P_BOLT_LMM_INF),]  
  stats <- subset(stats,stats$P_BOLT_LMM_INF<=5e-8)
  stats_1<-as.data.frame(stats)
  if(!is.na(stats_1[1,1])){
    myfiles_signal<-c(myfiles_signal,myfiles[i])
  }
}
  

len=length( myfiles_signal)
a=strsplit( myfiles_signal,'[/]')
filename=data.frame()
for (i in c(1:len)){
  c=a[[i]][9]   
  filename[i,1]=c
}




##choose principal components/traits with significant pvalue from mahattan plots 
SNP <- data.frame()
PC <- c()
for (i in 1:length(myfiles_signal)){
  stats <- fread(myfiles_signal[i],header = T)
  stats <- stats[order(stats$P_BOLT_LMM_INF),]  
  stats <- subset(stats,stats$P_BOLT_LMM_INF<=5e-8)
  for (j in 1:length(chr)){
    Locus <- subset(stats,stats$CHR == chr[j])
    SNP <- rbind(SNP,data.frame(c(Locus[1,c(1,2,3)])))  
    PC <- c(PC,filename[i,1])
  }
}
SNP <- cbind(SNP,PC)
colnames(SNP)[4] <- 'PC'   
write.table(file='sig_snps.txt',SNP,row.names = F,quote = F)

####prepare LD file for every significant snp  

for (i in 1:nrow(SNP)){
 #rsname-ÎªÎÄ¼þÃûÃüÃû
  cmd <- paste(plink_path,
               '--bfile',Uygur3M_path,
               '--r2 dprime',
               '--ld-snp',SNP$SNP[i],
               '--ld-window-r2',0,
               '--ld-window-kb',5000,
               '--ld-window',99999,
               '--out',paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/LD_file/',SNP$PC[i],'_',SNP$SNP[i],sep = '')
               )
  #--ld-window and --ld-window-r2 commands effectively means that output will be shown for all other SNPs within 1Mb of rs12345
  #adds the absolute value of Lewontin's D-prime statistic to table-formatted reports, and forces both r/r2 and D-prime to be based on the maximum likelihood solution to the cubic equation 
  #D=P(AB)-P(A)P(B)  r2=D * D / P(A)P(a)P(B)P(b)
  system(cmd) 
  ld <- read.table(paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/LD_file/',SNP$PC[i],'_',SNP$SNP[i],'.ld',sep = ''),header = T)
  ld_filter <- ld[,c(6,3,8,7)] #È¡½á¹ûµÄÌØ¶¨ÁÐ
  colnames(ld_filter) <- c('snp1','snp2','dprime','rsquare')
  write.table(file=paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/LD_file/',SNP$PC[i],'_',SNP$SNP[i],'.ld',sep = ''),ld_filter,row.names = F,quote = F)
}

### LocusZoom plot--demand on PYTHON27 
for (i in 1:nrow(SNP)){
  cmd <- paste(locuszoom_path,
               '--plotonly','--verbose',
               '--metal',paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/GWAS_result/bolt_assoc_resultprocess1/',SNP$PC[i],sep=''), #gwas½á¹û(without covariates)
               '--refsnp',SNP$SNP[i],
               '--chr',SNP$CHR[i],
               '--start',SNP$BP[i]-500000,
               '--end',SNP$BP[i] + 500000,
               '--prefix',paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/LocusZoom/',SNP$PC[i],'_',SNP$SNP[i],sep = ''), #output
               '--build hg19',
               '--ld',paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/LD_file/',SNP$PC[i],'_',SNP$SNP[i],'.ld',sep = ''),
               '--pvalcol P_BOLT_LMM_INF','--markercol SNP',
               'wightCol=Weight',
               '--snpset NULL','geneFontize=.8','smallDot=.3','largeDot=.9',
               'format=pdf','legened=auto','showRecomb=TRUE'
               )
  system(cmd)
}


## PAINTOR for fine mapping
if (!file.exists('FineMappingRunDirectory')){
  dir.create('FineMappingRunDirectory')
}
if (!file.exists('FineMappingOutDirectory')){
  dir.create('FineMappingOutDirectory')
}

## inpu.files
locus <- c()
for (i in 1:nrow(SNP)){
  locus <- c(locus,paste(SNP$PC[i],'_',SNP$SNP[i],sep = ''))
}
write.table(file='/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/FineMappingRunDirectory/inpu.files',locus,row.names = F,col.names = F,quote=F)

#SNP <- fread('sig_snps.txt',header = T)

####prepare LD file in matrix format and Zscore file 
for (i in 1:nrow(SNP)){
  stats <- fread(paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/GWAS_result/bolt_assoc_resultprocess1/',SNP$PC[i],sep=''),header = T)
  ZSCORE <- stats$BETA/stats$SE  #±ê×¼»¯µÄbeta
  PC <- data.frame(stats,ZSCORE)
  filter <- PC[PC$CHR==SNP$CHR[i]&PC$BP<=(SNP$BP[i]+500000)&PC$BP>=(SNP$BP[i]-500000),c(2,3,1,12)] #Ç°ºó50k,ÓëlocuszoomÑ¡¶¨¹Û²ì·¶Î§ÏàÍ¬
  write.table(file=paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/FineMappingRunDirectory/',SNP$PC[i],'_',SNP$SNP[i],sep = ''),filter,row.names = F,quote=F)
  write.table(file=paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/FineMappingRunDirectory/',SNP$PC[i],'_',SNP$SNP[i],'.snp',sep = ''),filter$SNP,col.names = F,row.names = F,quote=F)
  ## LD matrix file preparetion
  cmd <- paste(plink_path,
               '--bfile',Uygur3M_path,
               '--extract',paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/FineMappingRunDirectory/',SNP$PC[i],'_',SNP$SNP[i],'.snp',sep = ''), #½«Ñ¡¶¨µÄSNP set£¨Ç°ºó500kb) Á½Á½½øÐÐLDµÃµ½matrix
               '--r2 --matrix',
               '--out',paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/FineMappingRunDirectory/',SNP$PC[i],'_',SNP$SNP[i],sep = '')
  )
  system(cmd)
}

###prepare dummy annotation file
for (i in 1:nrow(SNP)){
  cmd <- paste('python2','/asnas/chenhua_group/guolu/toolkit/PAINTOR_V3.0-master/PAINTOR_Utilities/AnnotateLocus.py',
               '--input','/asnas/chenhua_group/guolu/2016_2018/dongyi_data/annotation_path.txt', #×¢ÊÍÎÄ¼þÂ·¾¶
               '--locus',paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/FineMappingRunDirectory/',SNP$PC[i],'_',SNP$SNP[i],sep = ''),
               '--out',paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/FineMappingRunDirectory/',SNP$PC[i],'_',SNP$SNP[i],'.annotations',sep = ''),
               '--chr CHR',
               '--pos BP')
  system(cmd)
  ## replace all value of 1 to 0 because we don't want any annotations
  annotations <- read.table(paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/FineMappingRunDirectory/',SNP$PC[i],'_',SNP$SNP[i],'.annotations',sep = ''),header = T)
  annotations[annotations$enhancers==1,1] = 0
  annotations[annotations$TSS==1,2] = 0
  write.table(file=paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/FineMappingRunDirectory/',SNP$PC[i],'_',SNP$SNP[i],'.annotations',sep = ''),annotations,row.names = F,quote = F)
}

## RUN PAINTOR
cmd <- paste(Paintor_path,
             '-input /asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/FineMappingRunDirectory/inpu.files',
             '-Zhead ZSCORE',
             '-LDname ld',
             '-in /asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/FineMappingRunDirectory/',
             '-out /asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/FineMappingOutDirectory/',
             '-mcmc',
             '-annotations enhancers,TSS')
system(cmd)
## change colname BP to pos
for (i in 1:nrow(SNP)){
  PAINTOR_result <- fread(paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/FineMappingOutDirectory/',SNP$PC[i],'_',SNP$SNP[i],'.results',sep = ''),header = T)
  colnames(PAINTOR_result)[2] <- 'pos'
  write.table(file=paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/FineMappingOutDirectory/',SNP$PC[i],'_',SNP$SNP[i],'.results',sep = ''),PAINTOR_result,row.names = F,quote = F)
}

## visualization with PAINTOR CANVIS
if (!file.exists('FineMappingResults')){
  dir.create('FineMappingResults')
}
for (i in 1:nrow(SNP)){
  cmd <- paste('python3 /asnas/chenhua_group/guolu/toolkit/PAINTOR_V3.0-master/CANVIS/CANVIS.py',
               '-l',paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/FineMappingOutDirectory/',SNP$PC[i],'_',SNP$SNP[i],'.results',sep = ''),
               '-z ZSCORE',
               '-r',paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/FineMappingRunDirectory/',SNP$PC[i],'_',SNP$SNP[i],'.ld',sep = ''),
               
               '-o',paste('/asnas/chenhua_group/guolu/2016_2018/dongyi_data/Causal_SNP/FineMappingResults/',SNP$PC[i],'_',SNP$SNP[i],sep = '')
  )
  system(cmd)
}

