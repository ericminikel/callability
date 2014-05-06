setwd("~/d/sci/050callab")
options(stringsAsFactors=FALSE)
require(car)
require(scatterplot3d)

cc = read.table("all_cc.txt",header=TRUE)
cp = read.table("all_cp.txt",header=TRUE)

gcolor = '#FFA824' # aureoline yellow
ecolor = '#0D4F8B' # indigo dye

# cdf is a vector giving the cdf for values 0 to 500
# the mean is sum over n of n*(CDF(n-1)-CDF(n))
gte0index = 6 # index of the column containing the "greater than or equal to 0" values
cp$meancov = rep(0.0,dim(cp)[1])
for (n in 1:500) {
  cp$meancov = cp$meancov + n*(cp[,n+gte0index-1]-cp[,n+gte0index])
}
subcov = cp[cp$sid=='65T_CR_1',c(1,2,3,4,5,507)]
m = lm(meancov ~ ilist + bamlist + as.factor(minbq) + as.factor(minmq), data=subcov)
summary(m)

dev.off()
barplot(subcov$meancov)
install.packages("scatterplot3d")


# compare WES vs. WGS mean coverage for all samples
all_ex = cp$minbq==0 & cp$minmq==0 & cp$ilist=="be.bed"
to_plot = rbind(cp$meancov[all_ex & cp$bamlist=="wgs.list"],
                cp$meancov[all_ex & cp$bamlist=="wes.list"])
colnames(to_plot) = cp$sid[all_ex & cp$bamlist=="wgs.list"]
to_plot
barplot(to_plot,beside=TRUE,col=c(gcolor,ecolor),border=NA,
        main='Mean coverage of exons, WES vs. WGS')
legend('bottomleft',c("WGS","WES"),col=c(gcolor,ecolor),pch=15)
# mean coverage of exons is categorically lower for WGS than WES

c65T = cp$minbq==0 & cp$minmq==0 & cp$sid == '65T_CR_1'
to_plot = rbind(cp$meancov[c65T & cp$bamlist=="wgs.list"],
                cp$meancov[c65T & cp$bamlist=="wes.list"])
colnames(to_plot) = cp$ilist[c65T & cp$bamlist=="wgs.list"]
barplot(to_plot,beside=TRUE,col=c(gcolor,ecolor),border=NA,
        main='Mean coverage of target sets in 65T, WES vs. WGS')
legend('bottomleft',c("WGS","WES"),col=c(gcolor,ecolor),pch=15)

# how does mean coverage drop off according to MQ?
subs = cp$minbq==0 & cp$sid == '65T_CR_1'
plot(cp$minmq[subs],cp$meancov[subs],pch='.',xlab='Min MQ',ylab='Mean coverage')
for (bamlist in unique(cp$bamlist)) {
  for (ilist in unique(cp$ilist)) {
    if (bamlist=="wgs.list") {k = gcolor} else if (bamlist=="wes.list") {k=ecolor}
    subs2 = subs & cp$ilist==ilist & cp$bamlist==bamlist
    points(cp$minmq[subs2], cp$meancov[subs2], type='b', lwd=3, pch=18, col=k)
  }
}
# conclusion: you take a big hit going from minimum MAPQ of 0 to 1
# but almost no hit from 1 to 10 or 20

# how about according to BQ?
subs = cp$minmq==0 & cp$sid == '65T_CR_1'
plot(cp$minbq[subs],cp$meancov[subs],pch='.',xlab='Min BQ',ylab='Mean coverage')
for (bamlist in unique(cp$bamlist)) {
  for (ilist in unique(cp$ilist)) {
    if (bamlist=="wgs.list") {k = gcolor} else if (bamlist=="wes.list") {k=ecolor}
    subs2 = subs & cp$ilist==ilist & cp$bamlist==bamlist
    points(cp$minbq[subs2], cp$meancov[subs2], type='b', lwd=3, pch=18, col=k)
  }
}
# conclusion: these genomes, but not exomes, have a lot of coverage contributed by
# bases with BQ < 10

# scatterplot matrix
pairs(~meancov + minmq + minbq, data=cp)
scatterplotMatrix(~meancov + minmq + minbq | bamlist , data=cp, col=c(ecolor,gcolor))
scatterplotMatrix(~meancov + minmq + minbq | ilist , data=cp )
# mean coverage appears to barely fall off according to minmq and minbq


# how about % at > 20x vs. min mq
subs = cp$minbq==0 & cp$sid == '65T_CR_1'
plot(cp$minmq[subs],cp$gte_20[subs],pch='.',xlab='Min MAPQ',ylab='% at 20x+')
for (bamlist in unique(cp$bamlist)) {
  for (ilist in unique(cp$ilist)) {
    if (bamlist=="wgs.list") {k = gcolor} else if (bamlist=="wes.list") {k=ecolor}
    subs2 = subs & cp$ilist==ilist & cp$bamlist==bamlist
    points(cp$minmq[subs2], cp$gte_20[subs2], type='b', lwd=3, pch=18, col=k)
  }
}
# conclusion: neither takes much of a hit with increasing min MAPQ
# but WGS just has much more uniform coverage, as virtually 100% of exonic
# bases are covered at 20x or better

# how about % at > 20x vs. min bq
subs = cp$minmq==0 & cp$sid == '65T_CR_1'
plot(cp$minbq[subs],cp$gte_20[subs],pch='.',xlab='Min BQ',ylab='% at 20x+')
for (bamlist in unique(cp$bamlist)) {
  for (ilist in unique(cp$ilist)) {
    if (bamlist=="wgs.list") {k = gcolor} else if (bamlist=="wes.list") {k=ecolor}
    subs2 = subs & cp$ilist==ilist & cp$bamlist==bamlist
    points(cp$minbq[subs2], cp$gte_20[subs2], type='b', lwd=3, pch=18, col=k)
  }
}
# conclusion: same as for min MQ

# plot the CDFs of WGS and WES for diff samples
subs = cp$minmq == 0 & cp$minbq == 0 & cp$ilist=="be.bed"
dmax=80 # maximum depth to bother plotting
plot(c(0,dmax),c(0,1),pch=NA,xlab='depth',ylab='CDF',xaxs='i',yaxs='i',
     main='Cumulative distribution function of coverage\nWES vs. WGS, 11 samples',
     sub='Minimum MAPQ = 0, minimum BQ = 0, all exonic bases')
for (sid in unique(cp$sid)) {
  for (bamlist in unique(cp$bamlist)) {
    subs2 = subs & cp$sid==sid & cp$bamlist==bamlist
    if (bamlist=="wgs.list") {k = gcolor} else if (bamlist=="wes.list") {k=ecolor}
    points(0:dmax,cp[subs2,gte0index:(gte0index+dmax)], type='l', lwd=3, col=k)
  }
}
legend('bottomleft',c("WES","WGS"),pch=NA,lwd=3,col=c(ecolor,gcolor),bty='n')
abline(v=20,col='red')
abline(v=10,col='red')
abline(h=.95,col='red')

# plot the CDFs of WGS and WES for diff samples
subs = cp$minmq == 0 & cp$minbq == 0 & cp$ilist=="bm.bed"
dmax=80 # maximum depth to bother plotting
plot(c(0,dmax),c(0,1),pch=NA,xlab='depth',ylab='CDF',xaxs='i',yaxs='i',
     main='Cumulative distribution function of coverage\nWES vs. WGS, 11 samples',
     sub='Minimum MAPQ = 0, minimum BQ = 0, muscle disease gene exons')
for (sid in unique(cp$sid)) {
  for (bamlist in unique(cp$bamlist)) {
    subs2 = subs & cp$sid==sid & cp$bamlist==bamlist
    if (bamlist=="wgs.list") {k = gcolor} else if (bamlist=="wes.list") {k=ecolor}
    points(0:dmax,cp[subs2,gte0index:(gte0index+dmax)], type='l', lwd=3, col=k)
  }
}
legend('bottomleft',c("WES","WGS"),pch=NA,lwd=3,col=c(ecolor,gcolor),bty='n')
abline(v=20,col='red')



# plot the CDFs of WGS and WES for diff samples
subs = cp$minmq == 1 & cp$minbq == 1 & cp$ilist=="be.bed"
dmax=80 # maximum depth to bother plotting
plot(c(0,dmax),c(0,1),pch=NA,xlab='depth',ylab='CDF',xaxs='i',yaxs='i',
     main='Cumulative distribution function of coverage\nWES vs. WGS, 11 samples',
     sub='Minimum MAPQ = 1, minimum BQ = 1, all exonic bases')
for (sid in unique(cp$sid)) {
  for (bamlist in unique(cp$bamlist)) {
    subs2 = subs & cp$sid==sid & cp$bamlist==bamlist
    if (bamlist=="wgs.list") {k = gcolor} else if (bamlist=="wes.list") {k=ecolor}
    points(0:dmax,cp[subs2,gte0index:(gte0index+dmax)], type='l', lwd=3, col=k)
  }
}
legend('bottomleft',c("WES","WGS"),pch=NA,lwd=3,col=c(ecolor,gcolor),bty='n')
abline(v=20,col='red')

# plot the CDFs of WGS and WES for diff samples - for exons, qual 20+
subs = cp$minmq == 20 & cp$minbq == 10 & cp$ilist=="be.bed"
dmax=80 # maximum depth to bother plotting
plot(c(0,dmax),c(0,1),pch=NA,xlab='depth',ylab='CDF',xaxs='i',yaxs='i',
     main='Cumulative distribution function of coverage\nWES vs. WGS, 11 samples',
     sub='Minimum MAPQ = 20, minimum BQ = 10, all exonic bases')
for (sid in unique(cp$sid)) {
  for (bamlist in unique(cp$bamlist)) {
    subs2 = subs & cp$sid==sid & cp$bamlist==bamlist
    if (bamlist=="wgs.list") {k = gcolor} else if (bamlist=="wes.list") {k=ecolor}
    points(0:dmax,cp[subs2,gte0index:(gte0index+dmax)], type='l', lwd=3, col=k)
  }
}
legend('bottomleft',c("WES","WGS"),pch=NA,lwd=3,col=c(ecolor,gcolor),bty='n')
abline(v=20,col='red')
abline(v=10,col='red')
abline(h=.95,col='red')


subs = cp$minmq == 1 & cp$minbq == 1 & cp$ilist=="bm.bed"
dmax=80 # maximum depth to bother plotting
plot(c(0,dmax),c(0,1),pch=NA,xlab='depth',ylab='CDF',xaxs='i',yaxs='i',
     main='Cumulative distribution function of coverage\nWES vs. WGS, 11 samples',
     sub='Minimum MAPQ = 1, minimum BQ = 1, muscle disease gene exons')
for (sid in unique(cp$sid)) {
  for (bamlist in unique(cp$bamlist)) {
    subs2 = subs & cp$sid==sid & cp$bamlist==bamlist
    if (bamlist=="wgs.list") {k = gcolor} else if (bamlist=="wes.list") {k=ecolor}
    points(0:dmax,cp[subs2,gte0index:(gte0index+dmax)], type='l', lwd=3, col=k)
  }
}
legend('bottomleft',c("WES","WGS"),pch=NA,lwd=3,col=c(ecolor,gcolor),bty='n')
abline(v=20,col='red')



subs = cp$minmq == 1 & cp$minbq == 1 & cp$ilist=="be10.bed"
dmax=80 # maximum depth to bother plotting
plot(c(0,dmax),c(0,1),pch=NA,xlab='depth',ylab='CDF',xaxs='i',yaxs='i',
     main='Cumulative distribution function of coverage\nWES vs. WGS, 11 samples',
     sub='Minimum MAPQ = 1, minimum BQ = 1, 10 bp buffer around exons')
for (sid in unique(cp$sid)) {
  for (bamlist in unique(cp$bamlist)) {
    subs2 = subs & cp$sid==sid & cp$bamlist==bamlist
    if (bamlist=="wgs.list") {k = gcolor} else if (bamlist=="wes.list") {k=ecolor}
    points(0:dmax,cp[subs2,gte0index:(gte0index+dmax)], type='l', lwd=3, col=k)
  }
}
legend('bottomleft',c("WES","WGS"),pch=NA,lwd=3,col=c(ecolor,gcolor),bty='n')
abline(v=20,col='red')



subs = cp$minmq == 1 & cp$minbq == 1 & cp$ilist=="bm10.bed"
dmax=80 # maximum depth to bother plotting
plot(c(0,dmax),c(0,1),pch=NA,xlab='depth',ylab='CDF',xaxs='i',yaxs='i',
     main='Cumulative distribution function of coverage\nWES vs. WGS, 11 samples',
     sub='Minimum MAPQ = 1, minimum BQ = 1, 10 bp buffer around muscle disease gene exons')
for (sid in unique(cp$sid)) {
  for (bamlist in unique(cp$bamlist)) {
    subs2 = subs & cp$sid==sid & cp$bamlist==bamlist
    if (bamlist=="wgs.list") {k = gcolor} else if (bamlist=="wes.list") {k=ecolor}
    points(0:dmax,cp[subs2,gte0index:(gte0index+dmax)], type='l', lwd=3, col=k)
  }
}
legend('bottomleft',c("WES","WGS"),pch=NA,lwd=3,col=c(ecolor,gcolor),bty='n')
abline(v=20,col='red')


ccmat = as.matrix(cc[cc$ilist=='be.bed' & 
                             cc$bamlist=='wgs.list' & 
                             cc$minbq==0 & 
                             cc$minmq==0,6:506])
k= colorRampPalette(c("white","darkgreen"))(100)

dev.off()
image(1:11,1:100,ccmat[,1:100],col=k,
      xlab='Number of samples covered at given depth',
      ylab='Depth',
      main='% of target covered in X samples at Y depth \
      Whole genome')
abline(h=20,col='red')
par(fig=c(.80,.99,.97,.99),new=TRUE)
par(mar=c(.1,.1,.1,.1))
plot(1:100,rep(1,100),col=k,type='h',lwd=5,ylim=c(0,1),cex=.1,
     yaxt='n',xaxt='n',ylab='',xaxs='i',yaxs='i',
     xlab='Percent of bases covered at given depth in given number of samples')
axis(side=1,at=(0:10)*10,labels=paste((0:10)*10,"%"),lwd.ticks=0,cex.axis=.6,padj=-4)


ccmat = as.matrix(cc[cc$ilist=='be.bed' & 
                       cc$bamlist=='wes.list' & 
                       cc$minbq==0 & 
                       cc$minmq==0,6:506])
k= colorRampPalette(c("white","darkgreen"))(100)
dev.off()
image(1:11,1:100,ccmat[,1:100],col=k,
      xlab='Number of samples covered at given depth',
      ylab='Depth',
      main='% of target covered in X samples at Y depth \
      Whole exome')
abline(h=20,col='red')
par(fig=c(.80,.99,.97,.99),new=TRUE)
par(mar=c(.1,.1,.1,.1))
plot(1:100,rep(1,100),col=k,type='h',lwd=5,ylim=c(0,1),cex=.1,
     yaxt='n',xaxt='n',ylab='',xaxs='i',yaxs='i',
     xlab='Percent of bases covered at given depth in given number of samples')
axis(side=1,at=(0:10)*10,labels=paste((0:10)*10,"%"),lwd.ticks=0,cex.axis=.6,padj=-4)



# 
# 
# # give more descriptive names
# cumulcov$ilist = gsub("be.bed","Broad exons",cumulcov$ilist)
# cumulcov$ilist = gsub("be10.bed","Broad exons 10bp buffer",cumulcov$ilist)
# cumulcov$ilist = gsub("bm.bed","Muscle exons",cumulcov$ilist)
# cumulcov$ilist = gsub("bm10.bed","Muscle exons 10bp buffer",cumulcov$ilist)
# 
# cumulcov$bamlist = gsub("wgs.list","WGS",cumulcov$bamlist)
# cumulcov$bamlist = gsub("wes.list","WES",cumulcov$bamlist)
# 
# cumulcov$minmq = paste("MAPQ >=",cumulcov$minmq)
# cumulcov$minbq = paste("BQ >=",cumulcov$minbq)
# 
# ?axis
# dev.off()
# 
# rcare = 1:80 # rcare = "range we care about"
# # i.e. i don't care how many bases are covered at >= 81x
# pct_cov = apply(ccmat,MARGIN=2,FUN=sum)/sum(ccmat[,1])
# pct_cov = pct_cov[-1] # lop off the gte_0 so the indices match the depth
# head(pct_cov)
# par(yaxs='i')
# plot(rcare,pct_cov[rcare],type='l',lwd=3,col='red',
#      xlab='Coverage',ylab='Percent of target',yaxt='n',
#      ylim=c(0,1))
# abline(h=c(0,1))
# axis(side=2,at=(0:10)/10,labels=paste((0:10)*10,"%",sep=""),cex.axis=.8)
# 
# ?par
# plot(0:10,0:10,xaxt='n')
# axis(side=1,at=(0:10)*10,labels=paste((0:10)*10,"%"),lwd.ticks=0,cex.axis=.4,mex=-2,padj=-1)


# now look at specific variants that were differentially called
wes = read.table('wes.table',header=TRUE,sep='\t')
wgs = read.table('wgs.table',header=TRUE,sep='\t')

fixcolnames = function(colnamevector) {
  newcolnames = tolower(gsub("\\.","_",colnamevector))
  return (newcolnames)
}

colnames(wes) = fixcolnames(colnames(wes))
colnames(wgs) = fixcolnames(colnames(wgs))

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

require(sqldf)

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
# 39

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
# 12

inboth = sqldf("
select   count(*)
from     wgs g, wes e
where    e.chrom = g.chrom
and      e.pos = g.pos
and      e.ref = g.ref
and      e.alt = g.alt
;")

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

barplot(as.matrix(rbind(wesonly,(inboth-exactgtsame)/2,exactgtsame,(inboth-exactgtsame)/2,wgsonly)),
        horiz=FALSE,beside=FALSE,border=NA,xaxt='n',yaxt='n',
        col=c(ecolor,'#999999','#000000','#999999',gcolor))
mtext(side=4,at=c(0,total/2,total),col=c(ecolor,'#000000',gcolor),las=2,line=-5,
      text=c('Called in WES only','11 exact genotypes identical','Called in WGS only'))

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

sum(wes$filter=='PASS')/dim(wes)[1] # .88
sum(wgs$filter=='PASS')/dim(wgs)[1] # .998

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

# select a few genotypes for followup
wes_only_pass$homref = paste(wes_only_pass$ref,'/',wes_only_pass$ref,sep='')
wgs_only_pass$homref = paste(wgs_only_pass$ref,'/',wgs_only_pass$ref,sep='')

# read in 

# create a sample-locus table for IGV_plotter
sltable = data.frame(locsamp=character(),locus=character(),sample=character())
# create a Google Spreadsheet for IGV validation 
gstable = data.frame(wgsigv=character(),wesigv=character(),sample=character(),
                     real=character(),comments=character())
pathprefix = "file:///Users/eric/d/sci/050callab/igv/"
tbls = list(WES=wes_only_pass,WGS=wgs_only_pass)
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


# now expand analysis to whole exome, but without looking at individual
# IGV screenshots

wes_e = read.table('wes_e.table',header=TRUE,sep='\t')
wgs_e = read.table('wgs_e.table',header=TRUE,sep='\t')

fixcolnames = function(colnamevector) {
  newcolnames = tolower(gsub("\\.","_",colnamevector))
  return (newcolnames)
}

colnames(wes_e) = fixcolnames(colnames(wes_e))
colnames(wgs_e) = fixcolnames(colnames(wgs_e))

# note that multi-allelics are already split into multiple lines
# by GATK VariantsToTable
# now convert variants to minimal representation
wes_e[,c("pos","ref","alt")] = minrep_vectorized(wes_e$pos,wes_e$ref,wes_e$alt)
wgs_e[,c("pos","ref","alt")] = minrep_vectorized(wgs_e$pos,wgs_e$ref,wgs_e$alt)

# contcatenate all genotypes into one string
# cat(paste("wes_e$",colnames(wes_e)[grep("gt",colnames(wes_e))],sep="",collapse=",\n"))
wes_e$gtcat = paste(wes_e$x15e_dd_1_gt,wes_e$x16e_md_1_gt,wes_e$x17e_pd_1_gt,wes_e$x23h_lm_1_gt,wes_e$x24h_cm_1_gt,wes_e$x25h_jm_1_gt,wes_e$x65t_cr_1_gt,wes_e$x66t_ng_1_gt,wes_e$x67t_sr_1_gt,wes_e$x68t_dr_1_gt,wes_e$x69t_gg_1_gt)
wes_e$gtcat = paste(wes_e$x15e_dd_1_gt,
                    wes_e$x16e_md_1_gt,
                    wes_e$x17e_pd_1_gt,
                    wes_e$x23h_lm_1_gt,
                    wes_e$x24h_cm_1_gt,
                    wes_e$x25h_jm_1_gt,
                    wes_e$x65t_cr_1_gt,
                    wes_e$x66t_ng_1_gt,
                    wes_e$x67t_sr_1_gt,
                    wes_e$x68t_dr_1_gt,
                    wes_e$x69t_gg_1_gt)
# cat(paste("wgs_e$",colnames(wes_e)[grep("gt",colnames(wes_e))],sep="",collapse=",\n"))
wgs_e$gtcat = paste(wgs_e$x15e_dd_1_gt,
                    wgs_e$x16e_md_1_gt,
                    wgs_e$x17e_pd_1_gt,
                    wgs_e$x23h_lm_1_gt,
                    wgs_e$x24h_cm_1_gt,
                    wgs_e$x25h_jm_1_gt,
                    wgs_e$x65t_cr_1_gt,
                    wgs_e$x66t_ng_1_gt,
                    wgs_e$x67t_sr_1_gt,
                    wgs_e$x68t_dr_1_gt,
                    wgs_e$x69t_gg_1_gt)

wes_eonly = sqldf("
select   count(*)
                  from     wes_e e
                  where    not exists (
                  select   null
                  from     wgs_e g
                  where    e.chrom = g.chrom
                  and      e.pos = g.pos
                  and      e.ref = g.ref
                  and      e.alt = g.alt
                  )
                  ;")
wes_eonly

wgs_eonly = sqldf("
                  select   count(*)
                  from     wgs_e g
                  where    not exists (
                  select   null
                  from     wes_e e
                  where    e.chrom = g.chrom
                  and      e.pos = g.pos
                  and      e.ref = g.ref
                  and      e.alt = g.alt
                  )
                  ;")
wgs_eonly

inboth_e = sqldf("
               select   count(*)
               from     wgs_e g, wes_e e
               where    e.chrom = g.chrom
               and      e.pos = g.pos
               and      e.ref = g.ref
               and      e.alt = g.alt
               ;")
inboth_e

exactgtsame_e = sqldf("
                    select   count(*)
                    from     wgs_e g, wes_e e
                    where    e.chrom = g.chrom
                    and      e.pos = g.pos
                    and      e.ref = g.ref
                    and      e.alt = g.alt
                    and      e.gtcat = g.gtcat
                    ;")
exactgtsame_e

total = wes_eonly+wgs_eonly+inboth_e

barplot(as.matrix(rbind(wes_eonly,(inboth_e-exactgtsame_e)/2,exactgtsame_e,(inboth_e-exactgtsame_e)/2,wgs_eonly)),
        horiz=FALSE,beside=FALSE,border=NA,xaxt='n',yaxt='n',
        col=c(ecolor,'#999999','#000000','#999999',gcolor))
mtext(side=4,at=c(0,total/2,total),col=c(ecolor,'#000000',gcolor),las=2,line=-5,
      text=c('Called in wes_e only','11 exact genotypes identical','Called in wgs_e only'))

wes_eonly_pass = sqldf("
select   *
from     wes_e e
                       where    e.filter == 'PASS'
                       and      not exists (
                       select   null
                       from     wgs_e g
                       where    e.chrom = g.chrom
                       and      e.pos = g.pos
                       and      e.ref = g.ref
                       and      e.alt = g.alt
                       )
                       ;")
dim(wes_eonly_pass)[1]/wes_eonly


wgs_eonly_pass = sqldf("
                       select   *
                       from     wgs_e g
                       where    g.filter == 'PASS'
                       and      not exists (
                       select   null
                       from     wes_e e
                       where    e.chrom = g.chrom
                       and      e.pos = g.pos
                       and      e.ref = g.ref
                       and      e.alt = g.alt
                       )
                       ;")
dim(wgs_eonly_pass)[1]/wgs_eonly

hist(wgs_eonly_pass$dp)
max(wgs_eonly_pass$dp)
hist(wes_eonly_pass$dp)

plot(density(wgs_e$dp))
points(density(wgs_eonly_pass$dp),type='l',col='red')

plot(density(wes_e$dp))
points(density(wes_eonly_pass$dp),type='l',col='red')
