# load required libraries
library(caTools)

# load data
## average gc in 1 kb sliding windows
gc<-read.table("soil-C-SFA/wsip/binning.d/H2-16.89/H2-16.89.gc.1000.tab",header=F)

## read depth from bedtools
f1<-read.table("soil-C-SFA/wsip/binning.d/H2-16.89/H2F1_16.genomecov.perbase.tab",header=F)
f1<-read.table("soil-C-SFA/wsip/binning.d/H2-16.89/H2F1_16.genomecov.perbase.tab",header=F)
f2<-read.table("soil-C-SFA/wsip/binning.d/H2-16.89/H2F2_16.genomecov.perbase.tab",header=F)
f3<-read.table("soil-C-SFA/wsip/binning.d/H2-16.89/H2F3_16.genomecov.perbase.tab",header=F)
f4<-read.table("soil-C-SFA/wsip/binning.d/H2-16.89/H2F4_16.genomecov.perbase.tab",header=F)
f5<-read.table("soil-C-SFA/wsip/binning.d/H2-16.89/H2F5_16.genomecov.perbase.tab",header=F)
f6<-read.table("soil-C-SFA/wsip/binning.d/H2-16.89/H2F6_16.genomecov.perbase.tab",header=F)
f7<-read.table("soil-C-SFA/wsip/binning.d/H2-16.89/H2F7_16.genomecov.perbase.tab",header=F)
f8<-read.table("soil-C-SFA/wsip/binning.d/H2-16.89/H2F8_16.genomecov.perbase.tab",header=F)
f9<-read.table("soil-C-SFA/wsip/binning.d/H2-16.89/H2F9_16.genomecov.perbase.tab",header=F)
all.f<-read.table("soil-C-SFA/wsip/binning.d/H2-16.89/HT2.genomecov.perbase.tab",header=F)

# calculate read depth in 1kb sliding windows
f1.1k.avg<-runmean(f1$V3,1000)
x1000<-seq(1,length(f1.1k.avg),1000)
f2.1k.avg<-runmean(f2$V3,1000)
f3.1k.avg<-runmean(f3$V3,1000)
f4.1k.avg<-runmean(f4$V3,1000)
f5.1k.avg<-runmean(f5$V3,1000)
f6.1k.avg<-runmean(f6$V3,1000)
f7.1k.avg<-runmean(f7$V3,1000)
f8.1k.avg<-runmean(f8$V3,1000)
f9.1k.avg<-runmean(f9$V3,1000)
all.f.1k.avg<-runmean(all.f$V3,1000)

# %GC of whole genome:
gc.total<-0.6794

# calculate average coverage for each library
f1.avg<-mean(f1$V3)
f2.avg<-mean(f2$V3)
f3.avg<-mean(f3$V3)
f4.avg<-mean(f4$V3)
f5.avg<-mean(f5$V3)
f6.avg<-mean(f6$V3)
f7.avg<-mean(f7$V3)
f8.avg<-mean(f8$V3)
f9.avg<-mean(f9$V3)
all.f.avg<-mean(all.f$V3)

png("soil-C-SFA/wsip/binning.d/H2-16.89/H2F1_16.genomecov.panel-graph.png",height=14,width=9,units="in",res=300,bg="transparent")
# frame dimensions
par(mfrow=c(12,1),mar=c(0.1,1,0.5,0.1)+0.1,
  oma=c(1,5,5,1)+0.1,
  xpd=NA)

# densest
plot(f1.1k.avg[x1000],pch=20,ylab="F1",xlab="",xaxt="n",cex.lab=2,cex.axis=1.5)
segments(x0=min(x1000),x1=(max(x1000)/1000),y0=f1.avg,y1=f1.avg,col="red",lwd=3)
#abline(h=f1.avg,col="red",lwd=3)

plot(f2.1k.avg[x1000],pch=20,ylab="F2",xlab="",xaxt="n",cex.lab=2,cex.axis=1.5)
segments(x0=min(x1000),x1=(max(x1000)/1000),y0=f2.avg,y1=f2.avg,col="red",lwd=3)
#abline(h=f2.avg,col="red",lwd=3)

plot(f3.1k.avg[x1000],pch=20,ylab="F3",xlab="",xaxt="n",cex.lab=2,cex.axis=1.5)
segments(x0=min(x1000),x1=(max(x1000)/1000),y0=f3.avg,y1=f3.avg,col="red",lwd=3)
#abline(h=f3.avg,col="red",lwd=3)

plot(f4.1k.avg[x1000],pch=20,ylab="F4",xlab="",xaxt="n",cex.lab=2,cex.axis=1.5)
segments(x0=min(x1000),x1=(max(x1000)/1000),y0=f4.avg,y1=f4.avg,col="red",lwd=3)
#abline(h=f4.avg,col="red",lwd=3)

plot(f5.1k.avg[x1000],pch=20,ylab="F5",xlab="",xaxt="n",cex.lab=2,cex.axis=1.5)
segments(x0=min(x1000),x1=(max(x1000)/1000),y0=f5.avg,y1=f5.avg,col="red",lwd=3)
#abline(h=f5.avg,col="red",lwd=3)

plot(f6.1k.avg[x1000],pch=20,ylab="F6",xlab="",xaxt="n",cex.lab=2,cex.axis=1.5)
segments(x0=min(x1000),x1=(max(x1000)/1000),y0=f6.avg,y1=f6.avg,col="red",lwd=3)
#abline(h=f6.avg,col="red",lwd=3)

plot(f8.1k.avg[x1000],pch=20,ylab="F8",xlab="",xaxt="n",cex.lab=2,cex.axis=1.5)
segments(x0=min(x1000),x1=(max(x1000)/1000),y0=f8.avg,y1=f8.avg,col="red",lwd=3)
#abline(h=f8.avg,col="red",lwd=3)

# least dense
plot(f9.1k.avg[x1000],pch=20,ylab="F9",xlab="",xaxt="n",cex.lab=2,cex.axis=1.5)
segments(x0=min(x1000),x1=(max(x1000)/1000),y0=f9.avg,y1=f9.avg,col="red",lwd=3)
#abline(h=f7.avg,col="red",lwd=3)

plot(all.f.1k.avg[x1000],pch=20,ylab="Unfractionated",xlab="",xaxt="n",cex.lab=2,cex.axis=1.5)
segments(x0=min(x1000),x1=(max(x1000)/1000),y0=all.f.avg,y1=all.f.avg,col="red",lwd=3)
#abline(h=all.f.avg,col="red",lwd=3)
#title(ylab="Coverage",line=5,outer=T)

#GC
plot(gc$V2,pch=20,ylab="%GC",xaxt="n",xlab="",cex.lab=2,cex.axis=1.5)
#title("Genome coordinates (1 kb)",line=+3)
title("Genome coordinates (1 kb)",line=-15,cex.main=2.5)
#axis(3,)
axis(1,cex.lab=2,cex.axis=1.5)
segments(x0=min(x1000),x1=(max(x1000)/1000),y0=gc.total,y1=gc.total,col="blue",lwd=3)
# abline(h=gc.total,col="blue",lwd=3)

#GC from scaffold one below
#plot(gc$V2[which(gc$V1=="H2-16_scaffold_1428274")],pch=20,ylab="%GC",xaxt="n",xlab="",cex.lab=2,cex.axis=1.5)
#title("Genome coordinates (1 kb)",line=-15,cex.main=2.5)
#axis(1,cex.lab=2,cex.axis=1.5)
#segments(x0=(min(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),x1=(max(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),y0=gc.total,y1=gc.total,col="blue",lwd=3)


dev.off()

png("soil-C-SFA/wsip/binning.d/H2-16.89/H2F1_16.genomecov.scaffold_1428274.panel-graph.png",height=15,width=9,units="in",res=300,bg="transparent")

### need to mess with this to have more blank space at the bottom. Currently it's at the top (where GC used to be)
par(mfrow=c(12,1),mar=c(0.1,1,0.5,0.1)+0.1,
  oma=c(1,5,5,1)+0.1,
  xpd=NA)

#densest
plot(f1.1k.avg[x1000[1666:1857]],pch=20,ylab="F1",xlab="",xaxt="n",cex.lab=2,cex.axis=1.5)
segments(x0=(min(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),x1=(max(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),y0=f1.avg,y1=f1.avg,col="red",lwd=3)

plot(f2.1k.avg[x1000[1666:1857]],pch=20,ylab="F2",xlab="",xaxt="n",cex.lab=2,cex.axis=1.5)
segments(x0=(min(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),x1=(max(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),y0=f2.avg,y1=f2.avg,col="red",lwd=3)

plot(f3.1k.avg[x1000[1666:1857]],pch=20,ylab="F3",xlab="",xaxt="n",cex.lab=2,cex.axis=1.5)
segments(x0=(min(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),x1=(max(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),y0=f3.avg,y1=f3.avg,col="red",lwd=3)

plot(f4.1k.avg[x1000[1666:1857]],pch=20,ylab="F4",xlab="",xaxt="n",cex.lab=2,cex.axis=1.5)
segments(x0=(min(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),x1=(max(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),y0=f4.avg,y1=f4.avg,col="red",lwd=3)

plot(f5.1k.avg[x1000[1666:1857]],pch=20,ylab="F5",xlab="",xaxt="n",cex.lab=2,cex.axis=1.5)
segments(x0=(min(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),x1=(max(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),y0=f5.avg,y1=f5.avg,col="red",lwd=3)

plot(f6.1k.avg[x1000[1666:1857]],pch=20,ylab="F6",xlab="",xaxt="n",cex.lab=2,cex.axis=1.5)
segments(x0=(min(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),x1=(max(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),y0=f6.avg,y1=f6.avg,col="red",lwd=3)

plot(f7.1k.avg[x1000[1666:1857]],pch=20,ylab="F7",xlab="",xaxt="n",cex.lab=2,cex.axis=1.5)
segments(x0=(min(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),x1=(max(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),y0=f7.avg,y1=f7.avg,col="red",lwd=3)

plot(f8.1k.avg[x1000[1666:1857]],pch=20,ylab="F8",xlab="",xaxt="n",cex.lab=2,cex.axis=1.5)
segments(x0=(min(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),x1=(max(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),y0=f8.avg,y1=f8.avg,col="red",lwd=3)

plot(f9.1k.avg[x1000[1666:1857]],pch=20,ylab="F9",xlab="",xaxt="n",cex.lab=2,cex.axis=1.5)
segments(x0=(min(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),x1=(max(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),y0=f9.avg,y1=f9.avg,col="red",lwd=3)

#unfragmented
plot(all.f.1k.avg[x1000[1666:1857]],pch=20,ylab="Unfractionated",xlab="",xaxt="n",cex.lab=2,cex.axis=1.5)
segments(x0=(min(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),x1=(max(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),y0=all.f.avg,y1=all.f.avg,col="red",lwd=3)

#GC
plot(gc$V2[which(gc$V1=="H2-16_scaffold_1428274")],pch=20,ylab="%GC",xaxt="n",xlab="",cex.lab=2,cex.axis=1.5)
title("Genome coordinates (1 kb)",line=-15,cex.main=2.5)
axis(1,cex.lab=2,cex.axis=1.5)
segments(x0=(min(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),x1=(max(gc$V3[which(gc$V1=="H2-16_scaffold_1428274")])/1000),y0=gc.total,y1=gc.total,col="blue",lwd=3)

dev.off()


#plot(gc$V2[which(gc$V1=="H2-16_scaffold_1428274")],pch=20,ylim=c(0.55,0.75),xlab="",ylab="%GC")
#plot(f1.1k.avg[x1000[1666:1857]],pch=20,ylab="F1 Cov",xlab="")
#plot(f2.1k.avg[x1000[1666:1857]],pch=20,ylab="F2 Cov",xlab="")
#plot(f3.1k.avg[x1000[1666:1857]],pch=20,ylab="F3 Cov",xlab="")
#plot(f4.1k.avg[x1000[1666:1857]],pch=20,ylab="F4 Cov",xlab="")
#plot(f5.1k.avg[x1000[1666:1857]],pch=20,ylab="F5 Cov",xlab="")
#plot(f6.1k.avg[x1000[1666:1857]],pch=20,ylab="F6 Cov",xlab="")
#plot(f7.1k.avg[x1000[1666:1857]],pch=20,ylab="F7 Cov",xlab="")
#plot(f8.1k.avg[x1000[1666:1857]],pch=20,ylab="F8 Cov",xlab="")
#plot(f9.1k.avg[x1000[1666:1857]],pch=20,ylab="F9 Cov",xlab="")
#plot(all.f.1k.avg[x1000[1666:1857]],pch=20,ylab="Unfractionated\n Cov",xlab="")
#plot(all.f.1k.avg[x1000[1666:1857]],pch=20,ylab="Unfractionated",xlab="")
#plot(all.f.1k.avg[x1000[1666:1857]],pch=20,ylab="Unfractionated",xlab="")
