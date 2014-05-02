setwd("~/d/sci/050callab")
options(stringsAsFactors=FALSE)

cumulcov = read.table("test1.sample_cumulative_coverage_counts")
ccmat = as.matrix(cumulcov)
k= colorRampPalette(c("white","darkgreen"))(100)

dev.off()
image(1:11,1:100,ccmat[,1:100],col=k,
      xlab='Number of samples covered at given depth',
      ylab='Depth',
      main='% of target covered in X samples at Y depth \
            Exome, MAPQ > 0, BQ > 20')
abline(h=20,col='black')
par(fig=c(.80,.99,.97,.99),new=TRUE)
par(mar=c(.1,.1,.1,.1))
plot(1:100,rep(1,100),col=k,type='h',lwd=5,ylim=c(0,1),cex=.1,
     yaxt='n',xaxt='n',ylab='',
     xlab='Percent of bases covered at given depth in given number of samples')
axis(side=1,at=(0:10)*10,labels=paste((0:10)*10,"%"),lwd.ticks=0,cex.axis=.6,padj=-4)

?axis
dev.off()

rcare = 1:80 # rcare = "range we care about"
# i.e. i don't care how many bases are covered at >= 81x
pct_cov = apply(ccmat,MARGIN=2,FUN=sum)/sum(ccmat[,1])
pct_cov = pct_cov[-1] # lop off the gte_0 so the indices match the depth
head(pct_cov)
par(yaxs='i')
plot(rcare,pct_cov[rcare],type='l',lwd=3,col='red',
     xlab='Coverage',ylab='Percent of target',yaxt='n',
     ylim=c(0,1))
abline(h=c(0,1))
axis(side=2,at=(0:10)/10,labels=paste((0:10)*10,"%",sep=""),cex.axis=.8)

?par
plot(0:10,0:10,xaxt='n')
axis(side=1,at=(0:10)*10,labels=paste((0:10)*10,"%"),lwd.ticks=0,cex.axis=.4,mex=-2,padj=-1)
