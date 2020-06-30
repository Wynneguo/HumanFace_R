#admixture plot
#input merge.Q_with_id | racefile.txt
#output MDS.pdf

data=read.table("C:/Users/GuoLu/Desktop/admix_merge.Q_with_id")
colnames(data)<-c('FID','IID','Q1','Q2')
#plot_1
barplot(t(as.matrix(data[,3:4])),space=c(0), col=rainbow(4),xlab="Individual #", ylab="Ancestry", border=NA) 
#space设置列间距
#plot_2 MDS
race<- read.table(file="C:/Users/GuoLu/Desktop/racefile.txt",header=TRUE)
datafile<- merge(data,race,by=c("IID","FID"))
head(datafile)

pdf("C:/Users/GuoLu/Desktop/MDS.pdf",width=7,height=7)
for (i in 1:nrow(datafile))
{
  if (datafile[i,5]=="EUR") {plot(datafile[i,3],datafile[i,4],type="p",xlim=c(0,1),ylim=c(0,1),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="green")
    par(new=T)}
  # par(new=T) 下一个函数不清空当前绘图
  if (datafile[i,5]=="ASN") {plot(datafile[i,3],datafile[i,4],type="p",xlim=c(0,1),ylim=c(0,1),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="red")
    par(new=T)}
  if (datafile[i,5]=="OWN") {plot(datafile[i,3],datafile[i,4],type="p",xlim=c(0,1),ylim=c(0,1),xlab="MDS Component 1",ylab="MDS Component 2",pch=3,cex=0.7,col="black") 
    par(new=T)}
}

abline(v=-0.035,lty=3)
abline(h=0.035,lty=3)
legend("topright", pch=c(1,1,1,1,3),c("EUR","ASN","AMR","AFR","OWN"),col=c("green","red",470,"blue","black"),bty="o",cex=1)
dev.off()
