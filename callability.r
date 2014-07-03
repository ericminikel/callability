setwd("~/d/sci/050callab")
options(stringsAsFactors=FALSE)
require(car)
require(sqldf)
require(VennDiagram)

cc = read.table("all_cc.txt",header=TRUE)
cp = read.table("all_cp.txt",header=TRUE)

gcolor = '#FFA824' # aureoline yellow
ecolor = '#0D4F8B' # indigo dye
# # average of the two colors
# bcolor = '#867C58'

# compute mean coverage.
# cdf is a vector giving the cdf for values 0 to 500
# the mean is sum over n of n*(CDF(n-1)-CDF(n))
gte0index = 6 # index of the column containing the "greater than or equal to 0" values
cp$meancov = rep(0.0,dim(cp)[1])
for (n in 1:500) {
  cp$meancov = cp$meancov + n*(cp[,n+gte0index-1]-cp[,n+gte0index])
}

# map sample names to fake sample names
name_mapping = data.frame(sid=unique(cp$sid),fakesid=paste("S",1:11,sep=""))
# map WGS and WES to their respective colors
k_mapping = data.frame(bamlist=c("wgs.hard.list","wes.hard.list"),k=c(gcolor,ecolor))

# compare WES vs. WGS mean coverage for all samples
png('img/meancov.png',width=600,height=400)
all_ex = cp$minbq==0 & cp$minmq==0 & cp$ilist=="be.bed"
to_plot = rbind(cp$meancov[all_ex & cp$bamlist=="wgs.hard.list"],
                cp$meancov[all_ex & cp$bamlist=="wes.hard.list"])
# colnames(to_plot) = cp$sid[all_ex & cp$bamlist=="wgs.hard.list"]
colnames(to_plot) = paste("S",1:11,sep="")
to_plot
barplot(to_plot,beside=TRUE,col=c(gcolor,ecolor),border=NA,
       main ="",cex.names=1,ylab='Mean depth')
#       main='Mean coverage of exons, WES vs. WGS',ylab='Mean depth')
title(expression("Mean coverage of Broad Exome, " * phantom("WGS") * " vs. " * phantom("WES")), col.main="black",cex.main=2)
title(expression(phantom("Mean coverage of Broad Exome, ") * "WGS" * phantom(" vs. ") * phantom("WES")), col.main=gcolor,font=2,cex.main=2)
title(expression(phantom("Mean coverage of Broad Exome, WGS vs. ") * "WES"), col.main=ecolor,font=2,cex.main=2)
# legend('bottomleft',c("WGS","WES"),col=c(gcolor,ecolor),pch=15,cex=.8)
abline(h=0,col='#000000')
dev.off() 
# mean coverage of exons is categorically lower for WGS than WES

# simple linear model of which factors influence depth
m = lm(meancov ~ ilist + bamlist + as.factor(minbq) + as.factor(minmq), data=cp)
summary(m)

png('img/meancov_bytargetset.png',width=600,height=400)
plot(NA,NA,xlim=c(1,4),ylim=c(1,120),xaxt='n',ylab='Mean depth',xlab='',
     main='Mean coverage of 11 samples by target set')
subs = cp$minbq==0 & cp$minmq==0
for (sid in unique(cp$sid)) {
  points(1:4,cp$meancov[subs & cp$sid==sid & cp$bamlist=='wgs.hard.list'],col=gcolor,type='b',pch=18,lwd=3)
  points(1:4,cp$meancov[subs & cp$sid==sid & cp$bamlist=='wes.hard.list'],col=ecolor,type='b',pch=18,lwd=3)
}
axis(side=1,at=1:4,labels=c("Broad exome","+-10bp","muscle exons","+-10bp"),lwd.ticks=0)
legend('bottomleft',c("WGS","WES"),col=c(gcolor,ecolor),pch=18,cex=.8)
dev.off()

# check monotonicity over minmq
for (minbq in unique(cp$minbq)) {
  for (bamlist in unique(cp$bamlist)) {
    for (ilist in unique(cp$ilist)) {
      for (sid in unique(cp$sid)) {
        subs2 = cp$minbq==minbq & cp$ilist==ilist & cp$bamlist==bamlist & cp$sid==sid
        if (!(all(cp$meancov[subs2]==cummin(cp$meancov[subs2])))) {
          print(paste(bamlist,ilist,minbq,sid))
        }
      }
    }
  }
}

# inspect non-monotonic instances
ilist = 'be.bed'
bamlist='wes.hard.list'
minbq = 1
sid = '25H_JM_1'
subs2 = cp$minbq==minbq & cp$ilist==ilist & cp$bamlist==bamlist & cp$sid==sid
cp$meancov[subs2]
tempview = cp[subs2,]
# the non-monotonicity is very slight ( <0.5x difference) and appears to owe to rounding error in the average.
# if you look at any one depth threshold, e.g. gte_15, gte_30, etc. they are all monotonic.
# the below code confirms that coverage at any given depth is monotonic.
mat = t(as.matrix(tempview[,6:506]))
which(mat[,2] > mat[,1])
which(mat[,3] > mat[,2])
which(mat[,4] > mat[,3])



# check monotonicity over minbq
for (minmq in unique(cp$minbq)) {
  for (bamlist in unique(cp$bamlist)) {
    for (ilist in unique(cp$ilist)) {
      for (sid in unique(cp$sid)) {
        subs2 = cp$minmq==minmq & cp$ilist==ilist & cp$bamlist==bamlist & cp$sid==sid
        if (!(all(cp$meancov[subs2]==cummin(cp$meancov[subs2])))) {
          print(paste(bamlist,ilist,minmq,sid))
        }
      }
    }
  }
}
# only one instance

# inspect non-monotonic instances
ilist = 'bm.bed'
bamlist='wes.hard.list'
minmq = 0
sid = '25H_JM_1'
subs2 = cp$minmq==minmq & cp$ilist==ilist & cp$bamlist==bamlist & cp$sid==sid
cp$meancov[subs2]
cp$minbq[subs2]
tempview = cp[subs2,]
# once again, this is just due to rounding error. below code shows that any one depth is montonic
mat = t(as.matrix(tempview[,6:506]))
which(mat[,2] > mat[,1])
which(mat[,3] > mat[,2])
which(mat[,4] > mat[,3])


# after the original 2014-05-06 run I *had* found depths where the coverage was non-monotonic
# over minbq or minmq. these examples are actually direct from GATK:
#$ cat cov_bm10.bed_wgs.hard.list_1_0.sample_cumulative_coverage_proportions | grep ^66T | cut -f62
#.76
#$cat cov_bm10.bed_wgs.hard.list_1_1.sample_cumulative_coverage_proportions | grep ^66T | cut -f62
#0.77



# how does mean coverage drop off according to MQ?
png('img/meancov_by_mapq_and_targetset.png',width=600,height=400)
subs = cp$minbq==0 & cp$sid == name_mapping$sid[name_mapping$fakesid=='S11']
plot(cp$minmq[subs],cp$meancov[subs],pch='.',xlab='Min MAPQ',ylab='Mean depth',
     xlim=c(0,23),main='Mean coverage by minimum MAPQ and target set in S11')
for (bamlist in unique(cp$bamlist)) {
  for (ilist in unique(cp$ilist)) {
    for (sid in unique(cp$sid)) {
      if (bamlist=="wgs.hard.list") {k = gcolor} else if (bamlist=="wes.hard.list") {k=ecolor}
      subs2 = subs & cp$ilist==ilist & cp$bamlist==bamlist & cp$sid==sid
      if (!(all(cp$meancov[subs2]==cummin(cp$meancov[subs2])))) {
        print(paste(bamlist,ilist,sid))
      }
      points(cp$minmq[subs2], cp$meancov[subs2], type='b', lwd=3, pch=18, col=k)
    }
  }
}
subs_lab = subs & cp$minmq == 20
text(cp$minmq[subs_lab], cp$meancov[subs_lab], labels=cp$ilist[subs_lab], pos=4,
     col=k_mapping$k[match(k_mapping$bamlist,cp$bamlist[subs_lab])])
legend('bottomleft',c("WGS","WES"),col=c(gcolor,ecolor),pch=18,cex=.8)
dev.off()
# conclusion: you take a big hit going from minimum MAPQ of 0 to 1
# but almost no hit from 1 to 10 or 20


m = lm(meancov ~ minmq*bamlist, data=subset(cp, minbq==0 & minmq %in% c(0,1)))
summary(m)
# on average you lose 2.6/88 = 3% of depth by requiring MAPQ>0 in WES
# and (2.6-1.35)/(88-25) = 2% in WGS


# how about according to BQ?
png('img/meancov_by_bq_and_targetset.png',width=600,height=400)
subs = cp$minmq==0 & cp$sid == name_mapping$sid[name_mapping$fakesid=='S11']
plot(cp$minbq[subs],cp$meancov[subs],pch='.',xlab='Min BQ',ylab='Mean depth',
     xlim=c(0,23),
     main='Mean coverage by minimum base quality and target set in S11')
for (sid in unique(cp$sid)) {
  for (bamlist in unique(cp$bamlist)) {
    for (ilist in unique(cp$ilist)) {
      if (bamlist=="wgs.hard.list") {k = gcolor} else if (bamlist=="wes.hard.list") {k=ecolor}
      subs2 = subs & cp$ilist==ilist & cp$bamlist==bamlist & cp$sid==sid
      points(cp$minbq[subs2], cp$meancov[subs2], type='b', lwd=3, pch=18, col=k)
    }
  }
}
subs_lab = subs & cp$minbq == 20
text(cp$minbq[subs_lab], cp$meancov[subs_lab], labels=cp$ilist[subs_lab], pos=4,
     col=k_mapping$k[match(k_mapping$bamlist,cp$bamlist[subs_lab])])
legend('bottomleft',c("WGS","WES"),col=c(gcolor,ecolor),pch=18,cex=.8)
dev.off()
# conclusion: these genomes, but not exomes, have a lot of coverage contributed by
# bases with BQ < 10

# scatterplot matrix
pairs(~meancov + minmq + minbq, data=cp)
scatterplotMatrix(~meancov + minmq + minbq | bamlist , data=cp, col=c(ecolor,gcolor))
scatterplotMatrix(~meancov + minmq + minbq | ilist , data=cp )
# mean coverage appears to barely fall off according to minmq and minbq


# how about % of Broad Exome at > 20x vs. min mq
png('img/be20plus_vs.minmq.png',width=600,height=400)
subs = cp$minbq==0 & cp$ilist=='be.bed'
plot(cp$minmq[subs],cp$gte_20[subs],pch='.',xlab='Min MAPQ',ylab='% at 20x+ depth',
     ylim=c(.75,1),yaxs='i',xaxt='n',yaxt='n',
     main='% of Broad Exome covered at 20x+
     by sample, MAPQ and WGS vs. WES')
axis(side=2,at=seq(.75,1,.05),labels=paste(seq(.75,1,.05)*100,"%",sep=""))
axis(side=1,at=c(0,1,10,20),labels=c(0,1,10,20))
for (sid in unique(cp$sid)) {
  for (bamlist in unique(cp$bamlist)) {
    
      if (bamlist=="wgs.hard.list") {k = gcolor} else if (bamlist=="wes.hard.list") {k=ecolor}
      subs2 = subs & cp$bamlist==bamlist & cp$sid==sid
      points(cp$minmq[subs2], cp$gte_20[subs2], type='b', lwd=3, pch=18, col=k)
    
  }
}
legend('bottomright',c("WGS","WES"),col=c(gcolor,ecolor),pch=18,cex=.8)
dev.off()
# conclusion: neither takes much of a hit with increasing min MAPQ
# but WGS just has much more uniform coverage, as virtually 100% of exonic
# bases are covered at 20x or better

# is this also true of 1x?

# how about % of Broad Exome at > 20x vs. min mq
png('img/be01plus_vs.minmq.png',width=600,height=400)
subs = cp$minbq==0 & cp$ilist=='be.bed'
plot(cp$minmq[subs],cp$gte_1[subs],pch='.',xlab='Min MAPQ',ylab='% at 1x+ depth',
     ylim=c(.75,1),yaxs='i',xaxt='n',yaxt='n',
     main='% of Broad Exome covered at 1x+
     by sample, MAPQ and WGS vs. WES')
axis(side=2,at=seq(.75,1,.05),labels=paste(seq(.75,1,.05)*100,"%",sep=""))
axis(side=1,at=c(0,1,10,20),labels=c(0,1,10,20))
for (sid in unique(cp$sid)) {
  for (bamlist in unique(cp$bamlist)) {
    
    if (bamlist=="wgs.hard.list") {k = gcolor} else if (bamlist=="wes.hard.list") {k=ecolor}
    subs2 = subs & cp$bamlist==bamlist & cp$sid==sid
    points(cp$minmq[subs2], cp$gte_1[subs2], type='b', lwd=3, pch=18, col=k)
    
  }
}
legend('bottomright',c("WGS","WES"),col=c(gcolor,ecolor),pch=18,cex=.8)
dev.off()


# how about % at > 20x vs. min bq
png('img/be20plus_vs.minbq.png',width=600,height=400)
subs = cp$minmq==0 & cp$ilist == 'be.bed'
plot(cp$minbq[subs],cp$gte_20[subs],pch='.',xlab='Min BQ',ylab='% at 20x+ depth',
     ylim=c(.75,1),yaxs='i',xaxt='n',yaxt='n',
     main='% of Broad Exome covered at 20x+
     by sample, BQ and WGS vs. WES')
axis(side=2,at=seq(.75,1,.05),labels=paste(seq(.75,1,.05)*100,"%",sep=""))
axis(side=1,at=c(0,1,10,20),labels=c(0,1,10,20))
for (bamlist in unique(cp$bamlist)) {
  for (sid in unique(cp$sid)) {
    if (bamlist=="wgs.hard.list") {k = gcolor} else if (bamlist=="wes.hard.list") {k=ecolor}
    subs2 = subs & cp$sid==sid & cp$bamlist==bamlist
    points(cp$minbq[subs2], cp$gte_20[subs2], type='b', lwd=3, pch=18, col=k)
  }
}
legend('bottomleft',c("WGS","WES"),col=c(gcolor,ecolor),pch=18,cex=.8)
dev.off()
# conclusion: same as for min MQ

# plot the CDFs of WGS and WES for diff samples
png('img/cdf_bysample_nomin.png',width=600,height=400)
subs = cp$minmq == 0 & cp$minbq == 0 & cp$ilist=="be.bed"
dmax=80 # maximum depth to bother plotting
plot(c(0,dmax),c(0,1),pch=NA,xlab='depth',xaxs='i',yaxs='i',
     ylab='Percentage of bases at this depth or better',yaxt='n',
     main='Cumulative distribution of coverage\nWES vs. WGS, 11 samples',
     sub='Minimum MAPQ = 0, minimum BQ = 0, all exonic bases')
axis(side=2,at=(0:10)/10,labels=paste((0:10)*10,"%",sep=""))
for (sid in unique(cp$sid)) {
  for (bamlist in unique(cp$bamlist)) {
    subs2 = subs & cp$sid==sid & cp$bamlist==bamlist
    if (bamlist=="wgs.hard.list") {k = gcolor} else if (bamlist=="wes.hard.list") {k=ecolor}
    points(0:dmax,cp[subs2,gte0index:(gte0index+dmax)], type='l', lwd=3, col=k)
  }
}
legend('bottomleft',c("WES","WGS"),pch=NA,lwd=3,col=c(ecolor,gcolor),bty='n')
abline(v=20,col='red')
dev.off()

# same but with MAPQ >= 1 and BQ >= 20

png('img/cdf_bysample_mapq1_bq20.png',width=600,height=400)
subs = cp$minmq == 1 & cp$minbq == 20 & cp$ilist=="be.bed"
dmax=80 # maximum depth to bother plotting
plot(c(0,dmax),c(0,1),pch=NA,xlab='depth',xaxs='i',yaxs='i',
     ylab='Percentage of bases at this depth or better',yaxt='n',
     main='Cumulative distribution of coverage\nWES vs. WGS, 11 samples',
     sub='Minimum MAPQ = 1, minimum BQ = 20, all exonic bases')
axis(side=2,at=(0:10)/10,labels=paste((0:10)*10,"%",sep=""))
for (sid in unique(cp$sid)) {
  for (bamlist in unique(cp$bamlist)) {
    subs2 = subs & cp$sid==sid & cp$bamlist==bamlist
    if (bamlist=="wgs.hard.list") {k = gcolor} else if (bamlist=="wes.hard.list") {k=ecolor}
    points(0:dmax,cp[subs2,gte0index:(gte0index+dmax)], type='l', lwd=3, col=k)
  }
}
legend('bottomleft',c("WES","WGS"),pch=NA,lwd=3,col=c(ecolor,gcolor),bty='n')
abline(v=10,col='red')
abline(v=20,col='red')
dev.off()

# check 10x+ coverage
wgs_10 = cp$gte_10[cp$minmq == 1 & cp$minbq == 20 & cp$ilist=="be.bed" & cp$bamlist=='wgs.hard.list']
wes_10 = cp$gte_10[cp$minmq == 1 & cp$minbq == 20 & cp$ilist=="be.bed" & cp$bamlist=='wes.hard.list']
mean(wgs_10)
mean(wes_10)

png('img/cdf_bysample_mapq1_bq20_bm.png',width=600,height=400)
subs = cp$minmq == 1 & cp$minbq == 20 & cp$ilist=="bm.bed"
dmax=80 # maximum depth to bother plotting
plot(c(0,dmax),c(0,1),pch=NA,xlab='depth',xaxs='i',yaxs='i',
     ylab='Percentage of bases at this depth or better',yaxt='n',
     main='Cumulative distribution of coverage\nWES vs. WGS, 11 samples',
     sub='Minimum MAPQ = 1, minimum BQ = 20, muscle disease exons')
axis(side=2,at=(0:10)/10,labels=paste((0:10)*10,"%",sep=""))
for (sid in unique(cp$sid)) {
  for (bamlist in unique(cp$bamlist)) {
    subs2 = subs & cp$sid==sid & cp$bamlist==bamlist
    if (bamlist=="wgs.hard.list") {k = gcolor} else if (bamlist=="wes.hard.list") {k=ecolor}
    points(0:dmax,cp[subs2,gte0index:(gte0index+dmax)], type='l', lwd=3, col=k)
  }
}
legend('bottomleft',c("WES","WGS"),pch=NA,lwd=3,col=c(ecolor,gcolor),bty='n')
abline(v=10,col='red')
abline(v=20,col='red')
dev.off()



png('img/cdf_bysample_mapq1_bq20_be10.png',width=600,height=400)
subs = cp$minmq == 1 & cp$minbq == 20 & cp$ilist=="be10.bed"
dmax=80 # maximum depth to bother plotting
plot(c(0,dmax),c(0,1),pch=NA,xlab='depth',xaxs='i',yaxs='i',
     ylab='Percentage of bases at this depth or better',yaxt='n',
     main='Cumulative distribution of coverage\nWES vs. WGS, 11 samples',
     sub='Minimum MAPQ = 1, minimum BQ = 20, 10bp buffer around Broad Exome')
axis(side=2,at=(0:10)/10,labels=paste((0:10)*10,"%",sep=""))
for (sid in unique(cp$sid)) {
  for (bamlist in unique(cp$bamlist)) {
    subs2 = subs & cp$sid==sid & cp$bamlist==bamlist
    if (bamlist=="wgs.hard.list") {k = gcolor} else if (bamlist=="wes.hard.list") {k=ecolor}
    points(0:dmax,cp[subs2,gte0index:(gte0index+dmax)], type='l', lwd=3, col=k)
  }
}
legend('bottomleft',c("WES","WGS"),pch=NA,lwd=3,col=c(ecolor,gcolor),bty='n')
abline(v=10,col='red')
abline(v=20,col='red')
dev.off()

# fancy 3-D plots where color saturation is percent of 
# bases in the exome at Y depth in X samples. 
# these plots are too hard to explain so I didn't end up
# using them in the blog post.
png('img/3dcov_be_wgs.png',width=400,height=400)
ccmat = as.matrix(cc[cc$ilist=='be.bed' & 
                             cc$bamlist=='wgs.hard.list' & 
                             cc$minbq==20 & 
                             cc$minmq==1,6:506])
#k= colorRampPalette(c("white","darkgreen"))(100)
k= colorRampPalette(c("white",gcolor))(100)
image(1:11,1:100,pmax(ccmat[,1:100]/ccmat[,1]-.8,0),col=k,
      xlab='Number of samples covered at given depth',
      ylab='Depth',
      main='% of target covered in X samples at Y depth in WGS',
      sub='Min MAPQ = 1, min BQ = 20, Broad Exome')
abline(h=20,col='red')
par(fig=c(.80,.99,.97,.99),new=TRUE)
par(mar=c(.1,.1,.1,.1))
plot(80:100,rep(1,21),col=k[(1:100)%%5==0],type='h',lwd=5,ylim=c(0,1),cex=.1,
     yaxt='n',xaxt='n',ylab='',xaxs='i',yaxs='i',
     xlab='Percent of bases covered at given depth in given number of samples')
axis(side=1,at=(8:10)*10,labels=paste((8:10)*10,"%"),lwd.ticks=0,cex.axis=.6,padj=-4)
dev.off()

png('img/3dcov_be_wes.png',width=400,height=400)
ccmat = as.matrix(cc[cc$ilist=='be.bed' & 
                       cc$bamlist=='wes.hard.list' & 
                       cc$minbq==0 & 
                       cc$minmq==0,6:506])
#k= colorRampPalette(c("white","darkgreen"))(100)
k= colorRampPalette(c("white",ecolor))(100)
image(1:11,1:100,pmax(ccmat[,1:100]/ccmat[,1]-.8,0),col=k,
      xlab='Number of samples covered at given depth',
      ylab='Depth',
      main='% of target covered in X samples at Y depth in WES')
abline(h=20,col='red')
par(fig=c(.80,.99,.97,.99),new=TRUE)
par(mar=c(.1,.1,.1,.1))
plot(80:100,rep(1,21),col=k[(1:100)%%5==0],type='h',lwd=5,ylim=c(0,1),cex=.1,
     yaxt='n',xaxt='n',ylab='',xaxs='i',yaxs='i',
     xlab='Percent of bases covered at given depth in given number of samples')
axis(side=1,at=(0:10)*10,labels=paste((0:10)*10,"%"),lwd.ticks=0,cex.axis=.6,padj=-4)
dev.off()
# didn't end up using these plots in the blog post as their meaning
# is a bit difficult to explain, and the colors are a bit hard to read.



##### variant call analysis

# now look at specific variants that were differentially called
wes = read.table('wes.table',header=TRUE,sep='\t')
wgs = read.table('wgs.table',header=TRUE,sep='\t')

fixcolnames = function(colnamevector) {
  newcolnames = tolower(gsub("\\.","_",colnamevector))
  return (newcolnames)
}

colnames(wes) = fixcolnames(colnames(wes))
colnames(wgs) = fixcolnames(colnames(wgs))

# see http://www.cureffi.org/2014/04/24/converting-genetic-variants-to-their-minimal-representation/
source("~/d/sci/src/minimal_representation/minimal-representation.r")

# note that multi-allelics are already split into multiple lines
# by GATK VariantsToTable
# now convert variants to minimal representation
wes[,c("pos","ref","alt")] = minrep_vectorized(wes$pos,wes$ref,wes$alt)
wgs[,c("pos","ref","alt")] = minrep_vectorized(wgs$pos,wgs$ref,wgs$alt)

# contcatenate all genotypes into one string
# cat(paste("wes$",colnames(wes)[grep("gt",colnames(wes))],sep="",collapse=",\n"))
wes$gtcat = paste(wes$x15e_dd_1_gt,wes$x16e_md_1_gt,wes$x17e_pd_1_gt,wes$x23h_lm_1_gt,wes$x24h_cm_1_gt,wes$x25h_jm_1_gt,wes$x65t_cr_1_gt,wes$x66t_ng_1_gt,wes$x67t_sr_1_gt,wes$x68t_dr_1_gt,wes$x69t_gg_1_gt)
wes$gtcat = paste(wes$x15e_dd_1_gt,
                  wes$x16e_md_1_gt,
                  wes$x17e_pd_1_gt,
                  wes$x23h_lm_1_gt,
                  wes$x24h_cm_1_gt,
                  wes$x25h_jm_1_gt,
                  wes$x65t_cr_1_gt,
                  wes$x66t_ng_1_gt,
                  wes$x67t_sr_1_gt,
                  wes$x68t_dr_1_gt,
                  wes$x69t_gg_1_gt)
# cat(paste("wgs$",colnames(wgs)[grep("gt",colnames(wgs))],sep="",collapse=",\n"))
wgs$gtcat = paste(wgs$x15e_dd_1_gt,
                  wgs$x16e_md_1_gt,
                  wgs$x17e_pd_1_gt,
                  wgs$x23h_lm_1_gt,
                  wgs$x24h_cm_1_gt,
                  wgs$x25h_jm_1_gt,
                  wgs$x65t_cr_1_gt,
                  wgs$x66t_ng_1_gt,
                  wgs$x67t_sr_1_gt,
                  wgs$x68t_dr_1_gt,
                  wgs$x69t_gg_1_gt)

# count variants that are shared vs. unqiue to WES or WGS

wesonly = sqldf("
select   count(*)
from     wes e
where    not exists (
      select   null
      from     wgs g
      where    e.chrom = g.chrom
      and      e.pos = g.pos
      and      e.ref = g.ref
      and      e.alt = g.alt
      )
;")


wgsonly = sqldf("
select   count(*)
from     wgs g
where    not exists (
      select   null
      from     wes e
      where    e.chrom = g.chrom
      and      e.pos = g.pos
      and      e.ref = g.ref
      and      e.alt = g.alt
)
;")


inboth = sqldf("
select   count(*)
from     wgs g, wes e
where    e.chrom = g.chrom
and      e.pos = g.pos
and      e.ref = g.ref
and      e.alt = g.alt
;")

# of variants in both sets, how many had exact GT same for all 11 samples?
exactgtsame = sqldf("
select   count(*)
from     wgs g, wes e
where    e.chrom = g.chrom
and      e.pos = g.pos
and      e.ref = g.ref
and      e.alt = g.alt
and      e.gtcat = g.gtcat
;")

total = wesonly+wgsonly+inboth

# proportion passing in WES and WGS
sum(wes$filter=='PASS')/dim(wes)[1] # .80
sum(wgs$filter=='PASS')/dim(wgs)[1] # .98

wes_only_pass = sqldf("
select   *
from     wes e
where    e.filter == 'PASS'
and      not exists (
      select   null
      from     wgs g
      where    e.chrom = g.chrom
      and      e.pos = g.pos
      and      e.ref = g.ref
      and      e.alt = g.alt
      )
;")
dim(wes_only_pass)[1]/wesonly


wgs_only_pass = sqldf("
select   *
from     wgs g
where    g.filter == 'PASS'
and      not exists (
      select   null
      from     wes e
      where    e.chrom = g.chrom
      and      e.pos = g.pos
      and      e.ref = g.ref
      and      e.alt = g.alt
      )
;")
dim(wgs_only_pass)[1]/wgsonly

inboth_passcount = sqldf("
select   count(*)
from     wgs g, wes e
where    e.chrom = g.chrom
and      e.pos = g.pos
and      e.ref = g.ref
and      e.alt = g.alt
and      (g.filter=='PASS' and e.filter=='PASS')
;") 
inboth_passcount
inboth_passcount/inboth

inboth_pass = sqldf("
select   *
from     wgs g, wes e
where    e.chrom = g.chrom
and      e.pos = g.pos
and      e.ref = g.ref
and      e.alt = g.alt
and      (g.filter=='PASS' and e.filter=='PASS')
;") 

wgsonly = as.integer(wgsonly)
wesonly = as.integer(wesonly)
inboth = as.integer(inboth)
wgsonly+inboth
wesonly+inboth
exactgtsame

png('img/alt_allele_venn.png',width=600,height=400)
draw.pairwise.venn(wesonly+inboth,wgsonly+inboth,inboth,
                   fill=c(ecolor,gcolor),lwd=0,alpha=c(.6,.8),
                   category=c("WES","WGS"),cex=.7,cat.fontfamily='sans',
                   euler.d=TRUE,scaled=TRUE,fontfamily='sans')
dev.off()

png('img/alt_allele_venn_passonly.png',width=600,height=400)
draw.pairwise.venn(dim(wes_only_pass)[1]+as.integer(inboth_passcount),
                   dim(wgs_only_pass)[1]+as.integer(inboth_passcount),
                   as.integer(inboth_passcount),
                   fill=c(ecolor,gcolor),lwd=0,alpha=c(.6,.8),
                   category=c("WES","WGS"),cex=.7,cat.fontfamily='sans',
                   euler.d=TRUE,scaled=TRUE,fontfamily='sans')
dev.off()
# didn't end up using this in the blog post because
# the VennDiagram package randomly switches the order of the circles
# on you and that makes it not look consistent from plot to plot.

dim(wes)[1]
dim(wgs)[1]

# # manually create a rectangular venn diagram
# barplot(height=c(1,1,1),width=c(wesonly,inboth,wgsonly),col=c(ecolor,bcolor,gcolor),
#         border=NA,space=c(0,0,0),yaxt='n',cex.main=.9,
#         main='Exonic variant alleles called in WGS vs. WES samples')
# atmarks = c(wesonly/2,wesonly+inboth/2,wesonly+inboth+wgsonly/2)
# for (i in 1:3) { # "bold" text by repeating it
#   mtext(side=1,at=atmarks,text=c("WES","both","WGS"),col=c(ecolor,bcolor,gcolor))
#   mtext(side=1,at=atmarks,text=c(wesonly,inboth,wgsonly),col=c(ecolor,bcolor,gcolor),padj=2)
# }
# # didn't end up using this because its rectangularness makes it not
# # not look enough like a venn diagram. darn.

# of the variants that were called in both sets but without exact match on all genotypes,
# does the allele count tend to be higher in WES or WGS calls?
acwgsmore = sqldf("
select   count(*)
from     wgs g, wes e
where    e.chrom = g.chrom
and      e.pos = g.pos
and      e.ref = g.ref
and      e.alt = g.alt
and      e.gtcat <> g.gtcat
and      g.ac > e.ac
;")
acwgsmore

acwesmore = sqldf("
select   count(*)
from     wgs g, wes e
where    e.chrom = g.chrom
and      e.pos = g.pos
and      e.ref = g.ref
and      e.alt = g.alt
and      e.gtcat <> g.gtcat
and      g.ac < e.ac
;")
acwesmore

# at what proportion of such sites was AC greater in WGS?
acwgsmore / (acwgsmore + acwesmore)

# which samples have the most discordant calls?
indiv_var_calls = data.frame(sid=character(),wescount=integer(),wgscount=integer())
wes_only_pass$homref = paste(wes_only_pass$ref,'/',wes_only_pass$ref,sep='')
wgs_only_pass$homref = paste(wgs_only_pass$ref,'/',wgs_only_pass$ref,sep='')
gtcols = colnames(wes_only_pass)[grep("_gt",colnames(wes_only_pass))]
adcols = which(colnames(wes_only_pass) %in% gtcols) + 1
gqcols = which(colnames(wes_only_pass) %in% gtcols) + 3
for (gtcol in gtcols) {
  wescount = sum(!(wes_only_pass[,gtcol] == wes_only_pass[,"homref"]) &
                   !(wes_only_pass[,gtcol] == './.'))
  wgscount = sum(!(wgs_only_pass[,gtcol] == wgs_only_pass[,"homref"]) &
                   !(wgs_only_pass[,gtcol] == './.'))
  sid = toupper(substr(gtcol,2,nchar(gtcol)-3))
  indiv_var_calls = rbind(indiv_var_calls,cbind(sid,wescount,wgscount)[1,])
}

# how many of the PASS variants in the WES have an alt genotype
# with GQ > 20, i.e. how many would be likely to have PASSed
# if not for joint-calling with other samples in whom the
# variant was more obviously polymorphic?

# create vector of PASS variant site indices with at least 1 alt genotype call with GQ > 20
wes_gq20_sites = integer(0)
for (site in 1:dim(wes_only_pass)[1]) {
  alt_called = which(!(wes_only_pass[site,gtcols] %in% c("./.",wes_only_pass[site,"homref"])))
  alt_gqs = wes_only_pass[site,gqcols[alt_called]]
  if (max(alt_gqs) > 20) {
    wes_gq20_sites = c(wes_gq20_sites,site)
  }
}
dim(wes_only_pass)[1]
length(wes_gq20_sites)
length(wes_gq20_sites)/dim(wes_only_pass)[1] # proportion of sites with at least 1 genotype call GQ > 20

wgs_gq20_sites = integer(0)
for (site in 1:dim(wgs_only_pass)[1]) {
  alt_called = which(!(wgs_only_pass[site,gtcols] %in% c("./.",wgs_only_pass[site,"homref"])))
  alt_gqs = wgs_only_pass[site,gqcols[alt_called]]
  if (max(alt_gqs) > 20) {
    wgs_gq20_sites = c(wgs_gq20_sites,site)
  }
}
dim(wgs_only_pass)[1]
length(wgs_gq20_sites)
length(wgs_gq20_sites)/dim(wgs_only_pass)[1] # proportion of sites with at least 1 genotype call GQ > 20

wgs_goodqual = length(wgs_gq20_sites)
wes_goodqual = length(wes_gq20_sites)
both_pass = as.integer(inboth_passcount)

png('img/alt_allele_venn_goodqual.png',width=600,height=400)
draw.pairwise.venn(
                   wgs_goodqual+both_pass,
                   wes_goodqual+both_pass,
                   both_pass,
                   fill=c(gcolor,ecolor),lwd=0,alpha=c(.6,.8),
                   category=c("WGS","WES"),cex=.7,cat.fontfamily='sans',
                   euler.d=TRUE,scaled=TRUE,fontfamily='sans')
dev.off()
# nevermind, just too hard to get these Venns to look good.

# select a random subset of 10 WES only and 10 WGS only for IGV followup
# use only high quality sites, i.e. PASS and at least 1 genotype GQ>20
# on my second attempt I discovered that some variants had very low GQ
# and only PASSed due to (presumably) higher quality in other samples
# with which they were joint-called
set.seed(1234)
to_val_wgs = wgs_only_pass[sample(wgs_gq20_sites,size=10,replace=FALSE),]
to_val_wes = wes_only_pass[sample(wes_gq20_sites,size=10,replace=FALSE),]

# create a sample-locus table for IGV_plotter
sltable = data.frame(locsamp=character(),locus=character(),sample=character())
# create a Google Spreadsheet for IGV validation 
gstable = data.frame(wgsigv=character(),wesigv=character(),sample=character(),
                     real=character(),comments=character())
pathprefix = "file:///Users/eric/d/sci/050callab/igv3/"
tbls = list(WES=to_val_wes,WGS=to_val_wgs)
# for all PASS variants called only in WES or only in WGS...
for (seqtype in names(tbls)) {
  tbl = tbls[[seqtype]]
  # get the names of the genotype columns
  gtcols = colnames(tbl)[grep("_gt",colnames(tbl))]
  # loop over the alleles (rows) in the table
  for (site in 1:dim(tbl)[1]) {
    # locus is chr:pos+50
    locus = paste(tbl$chrom[site],":",tbl$pos[site],"+50",sep="")
    # find which samples have a genotype that is neither homref nor ./.
    samples = gtcols[!(tbl[site,gtcols] %in% c("./.",tbl$homref[site]))]
    # grab the proper sample name by capitalizing and stripping x and _gt
    samples = toupper(substr(samples,2,nchar(samples)-3))
    # create a sample locus id for filenames
    locsamp = paste(samples,tbl[site,"chrom"],tbl[site,"pos"],tbl[site,"alt"],sep="_")
    # for each variant called only in WES or only in WGS, we want to see
    # a WES screenshot _and_ a WGS screenshot
    for (prefix in names(tbls)) {
      samples_pre = paste(prefix,samples,sep="_")
      locsamp_pre = paste(prefix,locsamp,sep="_")
      sltable = rbind(sltable,cbind(locsamp_pre,locus,samples_pre)[1,])
    }
    wgsigv = paste(pathprefix,"WGS_",locsamp,".png",sep="")
    wesigv = paste(pathprefix,"WES_",locsamp,".png",sep="")
    gstable = rbind(gstable,cbind(wgsigv,wesigv,samples,"","")[1,])
  }
}

write.table(sltable,"sample_locus.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(gstable,"gspread_validation.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

write.table(to_val_wgs,"to_val_wgs.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(to_val_wes,"to_val_wes.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")


# PDF of depth in shared vs. WGS-only or WES-only variants.
# didn't end up using this in blog post.
hist(wgs_eonly_pass$dp)
max(wgs_eonly_pass$dp)
hist(wes_eonly_pass$dp)

plot(density(wgs_e$dp))
points(density(wgs_eonly_pass$dp),type='l',col='red')

plot(density(wes_e$dp))
points(density(wes_eonly_pass$dp),type='l',col='red')


# looking up specific sites from one dataset to see whether/why 
# they weren't called in the other dataset

wgs[wgs$chrom=='6' & wgs$pos %in% 31323864:31325864,1:9] # 31324864
wgs[wgs$chrom=='1' & wgs$pos == 152186222,]

dim(wgs[wgs$chrom=='10' & wgs$pos > 135491000 & wgs$pos < 135492000,])

dim(wes[wes$chrom=='10' & wes$pos > 135491000 & wes$pos < 135492000,])

# attempt to find this one variant at 1:171162540 in wgs calls
wgs[wgs$chrom=='1' & wgs$pos > 171162500 & wgs$pos < 171162600,]
# it is really not there

