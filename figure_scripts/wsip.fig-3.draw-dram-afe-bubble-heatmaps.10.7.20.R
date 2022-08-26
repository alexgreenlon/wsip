library(cowplot)
library(ggplot2)
library(gridExtra)
library(grid)
library(readxl)
library(svglite)

#dat<-read.csv("wsip/results/dram/metawrap-drep-bins/genome-summaries/metabolism_summary-plus-pathways.afe.median.sitewise-ks-tests.long.csv",sep='\t',header=T)
#dat<-read.csv("wsip/results/ks-stats/combined-annotations.fixed-afe-median.sitewise-ks-tests.fixed.csv",sep=',',header=T)
#dat<-read.csv("soil-C-SFA/wsip/results/dram/metawrap-drep-bins/genome-summaries/metabolism_summary.fixed-afe-median.sitewise-ks-tests.amended.tsv",sep='\t',header=T)
dat<-read.csv("wsip/results/ks-stats/combined-annotations.fixed-afe-median.sitewise-ks-tests.fixed.tsv",sep='\t',header=T)
dat<-as.data.frame(dat)

#bin.info<-read.csv("soil-C-SFA/wsip/binning.d/metawrap.bin_refinement_comp50_cont_25.4.2.20/metawrap-bins.ggkbase-organism_info.drep-set.txt",sep='\t',header=T)
bin.info<-read.csv("soil-C-SFA/docs/wsip-metawrap-bins.genome-info.10.7.20.csv",sep='\t',header=T)

####### normalize data so AFE values are scaled by range at each sites and bin counts are scaled by total bins at each site...

dat$A.afe.median.mean.norm<-dat$A.afe.median.mean/max(dat$A.afe.median.mean,na.rm=T)
dat$H.afe.median.mean.norm<-dat$H.afe.median.mean/max(dat$H.afe.median.mean,na.rm=T)
dat$S.afe.median.mean.norm<-dat$S.afe.median.mean/max(dat$S.afe.median.mean,na.rm=T)

dat$A.hits.norm<-dat$A.hits/max(dat$A.hits,na.rm=T)
dat$H.hits.norm<-dat$H.hits/max(dat$H.hits,na.rm=T)
dat$S.hits.norm<-dat$S.hits/max(dat$S.hits,na.rm=T)

#bin.info$site<-substr(bin.info$name,0,1)

site.afes<-aggregate(afe.median ~ site, data=bin.info[which(bin.info$drep.set==1),],FUN=mean)
afe.A<-site.afes[1,2]
afe.H<-site.afes[2,2]
afe.S<-site.afes[3,2]

site.min.afes<-aggregate(afe.median ~ site, data=bin.info[which(bin.info$drep.set==1),],FUN=min)
afe.min.A<-site.min.afes[1,2]
afe.min.H<-site.min.afes[2,2]
afe.min.S<-site.min.afes[3,2]

site.max.afes<-aggregate(afe.median ~ site, data=bin.info[which(bin.info$drep.set==1),],FUN=max)
afe.max.A<-site.max.afes[1,2]
afe.max.H<-site.max.afes[2,2]
afe.max.S<-site.max.afes[3,2]

afe.mean.A.norm<-afe.A/afe.max.A
afe.min.A.norm<-afe.min.A/afe.max.A
afe.mean.H.norm<-afe.H/afe.max.H
afe.min.H.norm<-afe.min.H/afe.max.H
afe.mean.S.norm<-afe.S/afe.max.S
afe.min.S.norm<-afe.min.S/afe.max.S

dat$A.sig.ks<-as.factor(ifelse(dat$A.p.adj<0.05,1,0))
dat$H.sig.ks<-as.factor(ifelse(dat$H.p.adj<0.05,1,0))
dat$S.sig.ks<-as.factor(ifelse(dat$S.p.adj<0.05,1,0))

### Consolidated (without EPS)

# then just need to fix the legend names (or could do in illustrator)
# and make one for each site

afe.mid<-mean(c(afe.mean.A.norm,afe.mean.H.norm,afe.mean.S.norm))
afe.hi<-max(c(dat[which(dat$bubble.plot.order>0),]$H.afe.median.mean.norm,dat[which(dat$bubble.plot.order>0),]$A.afe.median.mean.norm,dat[which(dat$bubble.plot.order>0),]$S.afe.median.mean.norm),na.rm=T)
afe.low<-min(c(dat[which(dat$bubble.plot.order>0),]$H.afe.median.mean.norm,dat[which(dat$bubble.plot.order>0),]$A.afe.median.mean.norm,dat[which(dat$bubble.plot.order>0),]$S.afe.median.mean.norm),na.rm=T)

bin.hi<-max(c(dat[which(dat$bubble.plot.order>0),]$A.hits.norm,dat[which(dat$bubble.plot.order>0),]$H.hits.norm,dat[which(dat$bubble.plot.order>0),]$S.hits.norm),na.rm=T)
bin.low<-min(c(dat[which(dat$bubble.plot.order>0),]$A.hits.norm,dat[which(dat$bubble.plot.order>0),]$H.hits.norm,dat[which(dat$bubble.plot.order>0),]$S.hits.norm),na.rm=T)


dat$bubble.plot.desc<-factor(dat$bubble.plot.desc, levels=dat$bubble.plot.desc[order(-dat$bubble.plot.order)])

p1<-ggplot(data=dat[which(dat$bubble.plot.order>0),],aes(x="Angelo",y=bubble.plot.desc)) +
  geom_point(aes(fill=A.afe.median.mean,size=A.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.A,limits=c(afe.min.A,afe.max.A)) +
  geom_point(data=dat[which(dat$bubble.plot.order>0),],aes(x="Angelo",y=bubble.plot.desc,color=A.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/combined-bubble-plots.angelo.pdf",width=5,height=5)

p1<-ggplot(data=dat[which(dat$bubble.plot.order>0),],aes(x="Hopland",y=bubble.plot.desc)) +
  geom_point(aes(fill=H.afe.median.mean,size=H.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.H,limits=c(afe.min.H,afe.max.H)) +
  geom_point(data=dat[which(dat$bubble.plot.order>0),],aes(x="Hopland",y=bubble.plot.desc,color=H.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/combined-bubble-plots.hopland.pdf",width=5,height=5)

p1<-ggplot(data=dat[which(dat$bubble.plot.order>0),],aes(x="Sedgwick",y=bubble.plot.desc)) +
  geom_point(aes(fill=S.afe.median.mean,size=S.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.S,limits=c(afe.min.S,afe.max.S)) +
  geom_point(data=dat[which(dat$bubble.plot.order>0),],aes(x="Sedgwick",y=bubble.plot.desc,color=S.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/combined-bubble-plots.sedgwick.pdf",width=5,height=5)

#### normalized

p1<-ggplot(data=dat[which(dat$bubble.plot.order>0),],aes(x="Angelo",y=bubble.plot.desc)) +
  geom_point(aes(fill=A.afe.median.mean.norm,size=A.hits.norm),pch=21,color="black") + scale_size(limits=c(bin.low,bin.hi)) +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0.5,limits=c(afe.low,afe.hi)) +
  geom_point(data=dat[which(dat$bubble.plot.order>0),],aes(x="Angelo",y=bubble.plot.desc,color=A.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/combined-bubble-plots.angelo.norm.png",dpi=300,width=11,height=10,units = "in")
ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/combined-bubble-plots.angelo.norm.svg",dpi=300,width=13.5,height=10,units = "in")

#### This is how I actually saved it, as of 9/6/21
pdf("wsip/results/ks-stats/pathway-bubble-plots/combined-bubble-plots.angelo.norm.pdf")
p1
dev.off()
# did this in rstudio with whatever the graphics window was set up as
# then converted it to png from preview to get into google slides :.-(

# otherwise playing around with the point and text sizes in ggplot, then using ggsave to export a png
# I run into too many issues with the various elements changing size as you scale the png, 
# or the resolution is too low

p2<-ggplot(data=dat[which(dat$bubble.plot.order>0),],aes(x="Hopland",y=bubble.plot.desc)) +
  geom_point(aes(fill=H.afe.median.mean.norm,size=H.hits.norm),pch=21,color="black") + scale_size(limits=c(bin.low,bin.hi)) +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0.5,limits=c(afe.low,afe.hi)) +
  geom_point(data=dat[which(dat$bubble.plot.order>0),],aes(x="Hopland",y=bubble.plot.desc,color=H.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

pdf("wsip/results/ks-stats/pathway-bubble-plots/combined-bubble-plots.hopland.norm.pdf")
p2
dev.off()


ggsave(p2,file="wsip/results/ks-stats/pathway-bubble-plots/combined-bubble-plots.hopland.norm.pdf",width=5,height=5)

p3<-ggplot(data=dat[which(dat$bubble.plot.order>0),],aes(x="Sedgwick",y=bubble.plot.desc)) +
  geom_point(aes(fill=S.afe.median.mean.norm,size=S.hits.norm),pch=21,color="black") + scale_size(limits=c(bin.low,bin.hi)) +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0.5,limits=c(afe.low,afe.hi)) +
  geom_point(data=dat[which(dat$bubble.plot.order>0),],aes(x="Sedgwick",y=bubble.plot.desc,color=S.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

pdf("wsip/results/ks-stats/pathway-bubble-plots/combined-bubble-plots.sedgwick.norm.pdf")
p3
dev.off()

ggsave(p3,file="wsip/results/ks-stats/pathway-bubble-plots/combined-bubble-plots.sedgwick.norm.pdf",width=5,height=5)

# Cutinase
# Angelo
p1<-ggplot(dat[which(dat$annotation.source=="DRAM" & dat$gene == "CE5"),],aes(x="Angelo",y=gene)) +
  geom_point(aes(fill=A.afe.median.mean,size=A.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.A,limits=c(afe.min.A,afe.max.A)) +
#  geom_point(data=dat[which(dat$annotation.source=="DRAM" & dat$gene == "CE5"),],aes(x="Angelo",y=gene,color=A.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
#  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/A.cutinase.png",width=5,height=7)

# Hopland
p1<-ggplot(dat[which(dat$annotation.source=="DRAM" & dat$gene == "CE5"),],aes(x="Hopland",y=gene)) +
  geom_point(aes(fill=H.afe.median.mean,size=H.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.H,limits=c(afe.min.H,afe.max.H)) +
#  geom_point(data=dat[which(dat$annotation.source=="DRAM" & dat$gene == "CE5"),],aes(x="Hopland",y=gene,color=H.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
#  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/H.cutinase.png",width=5,height=7)

# Sedgwick
p1<-ggplot(dat[which(dat$annotation.source=="DRAM" & dat$gene == "CE5"),],aes(x="Sedgwick",y=gene)) +
  geom_point(aes(fill=S.afe.median.mean,size=S.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.S,limits=c(afe.min.S,afe.max.S)) +
  geom_point(data=dat[which(dat$annotation.source=="DRAM" & dat$gene == "CE5"),],aes(x="Sedgwick",y=gene,color=S.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/S.cutinase.png",width=5,height=7)

# Flagellin
p1<-ggplot(dat[which(dat$annotation.source=="DRAM" & dat$gene == "K02406"),],aes(x="Angelo",y=gene)) +
  geom_point(aes(fill=A.afe.median.mean,size=A.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.A,limits=c(afe.min.A,afe.max.A)) +
  geom_point(data=dat[which(dat$annotation.source=="DRAM" & dat$gene == "K02406"),],aes(x="Angelo",y=gene,color=A.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/A.flagellin.png",width=5,height=7)

# Hopland
p1<-ggplot(dat[which(dat$annotation.source=="DRAM" & dat$gene == "K02406"),],aes(x="Hopland",y=gene)) +
  geom_point(aes(fill=H.afe.median.mean,size=H.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.H,limits=c(afe.min.H,afe.max.H)) +
  geom_point(data=dat[which(dat$annotation.source=="DRAM" & dat$gene == "K02406"),],aes(x="Hopland",y=gene,color=H.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/H.flagellin.png",width=5,height=7)

# Sedgwick
p1<-ggplot(dat[which(dat$annotation.source=="DRAM" & dat$gene == "K02406"),],aes(x="Sedgwick",y=gene)) +
  geom_point(aes(fill=S.afe.median.mean,size=S.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.S,limits=c(afe.min.S,afe.max.S)) +
#  geom_point(data=dat[which(dat$annotation.source=="DRAM" & dat$gene == "K02406"),],aes(x="Sedgwick",y=gene,color=S.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
#  scale_color_manual(values=c("black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/S.flagellin.png",width=5,height=7)

# DRAM CAZY pathways
# Angelo
p1<-ggplot(dat[which(dat$annotation.source=="DRAM product" & dat$header == "CAZY pathways"),],aes(x="Angelo",y=gene)) +
  geom_point(aes(fill=A.afe.median.mean,size=A.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.A,limits=c(afe.min.A,afe.max.A)) +
  geom_point(data=dat[which(dat$annotation.source=="DRAM product" & dat$header == "CAZY pathways"),],aes(x="Angelo",y=gene,color=A.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/A.cazy-paths.png",width=5,height=7)

# Hopland
p1<-ggplot(dat[which(dat$annotation.source=="DRAM product" & dat$header == "CAZY pathways"),],aes(x="Hopland",y=gene)) +
  geom_point(aes(fill=H.afe.median.mean,size=H.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.H,limits=c(afe.min.H,afe.max.H)) +
  geom_point(data=dat[which(dat$annotation.source=="DRAM product" & dat$header == "CAZY pathways"),],aes(x="Hopland",y=gene,color=H.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/H.cazy-paths.png",width=5,height=7)

# Sedgwick
p1<-ggplot(dat[which(dat$annotation.source=="DRAM product" & dat$header == "CAZY pathways"),],aes(x="Sedgwick",y=gene)) +
  geom_point(aes(fill=S.afe.median.mean,size=S.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.S,limits=c(afe.min.S,afe.max.S)) +
  geom_point(data=dat[which(dat$annotation.source=="DRAM product" & dat$header == "CAZY pathways"),],aes(x="Sedgwick",y=gene,color=S.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/S.cazy-paths.png",width=5,height=7)

# DRAM Nitrogen pathways
# Angelo
p1<-ggplot(dat[which(dat$annotation.source=="DRAM product" & dat$header == "Nitrogen.metabolism"),],aes(x="Angelo",y=gene)) +
  geom_point(aes(fill=A.afe.median.mean,size=A.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.A,limits=c(afe.min.A,afe.max.A)) +
  geom_point(data=dat[which(dat$annotation.source=="DRAM product" & dat$header == "Nitrogen.metabolism"),],aes(x="Angelo",y=gene,color=A.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/A.n-paths.png",width=7,height=7)

# Hopland
p1<-ggplot(dat[which(dat$annotation.source=="DRAM product" & dat$header == "Nitrogen.metabolism"),],aes(x="Hopland",y=gene)) +
  geom_point(aes(fill=H.afe.median.mean,size=H.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.H,limits=c(afe.min.H,afe.max.H)) +
  geom_point(data=dat[which(dat$annotation.source=="DRAM product" & dat$header == "Nitrogen.metabolism"),],aes(x="Hopland",y=gene,color=H.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/H.n-paths.png",width=7,height=7)

# Sedgwick
p1<-ggplot(dat[which(dat$annotation.source=="DRAM product" & dat$header == "Nitrogen.metabolism"),],aes(x="Sedgwick",y=gene)) +
  geom_point(aes(fill=S.afe.median.mean,size=S.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.S,limits=c(afe.min.S,afe.max.S)) +
  geom_point(data=dat[which(dat$annotation.source=="DRAM product" & dat$header == "Nitrogen.metabolism"),],aes(x="Sedgwick",y=gene,color=S.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/S.n-paths.png",width=7,height=7)

# EPS
# Angelo
p1<-ggplot(dat[which(dat$annotation.source=="EPS" & dat$header == "Synthase-dependent EPS"),],aes(x="Angelo",y=gene)) +
  geom_point(aes(fill=A.afe.median.mean,size=A.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.A,limits=c(afe.min.A,afe.max.A)) +
  geom_point(data=dat[which(dat$annotation.source=="EPS" & dat$header == "Synthase-dependent EPS"),],aes(x="Angelo",y=gene,color=A.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/A.eps.png",width=3.5,height=2.5)

# Hopland
p1<-ggplot(dat[which(dat$annotation.source=="EPS" & dat$header == "Synthase-dependent EPS"),],aes(x="Hopland",y=gene)) +
  geom_point(aes(fill=H.afe.median.mean,size=H.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.H,limits=c(afe.min.H,afe.max.H)) +
  geom_point(data=dat[which(dat$annotation.source=="EPS" & dat$header == "Synthase-dependent EPS"),],aes(x="Hopland",y=gene,color=H.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/H.eps.png",width=3.5,height=2.5)

# Sedgwick
p1<-ggplot(dat[which(dat$annotation.source=="EPS" & dat$header == "Synthase-dependent EPS"),],aes(x="Sedgwick",y=gene)) +
  geom_point(aes(fill=S.afe.median.mean,size=S.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.S,limits=c(afe.min.S,afe.max.S)) +
  geom_point(data=dat[which(dat$annotation.source=="EPS" & dat$header == "Synthase-dependent EPS"),],aes(x="Sedgwick",y=gene,color=S.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/S.eps.png",width=3.5,height=2.5)

# METABLIC c1 metabolism
# Angelo
p1<-ggplot(dat[which(dat$annotation.source=="METABOLIC" & dat$header == "C1 metabolism"),],aes(x="Angelo",y=paste(gene_description,module,sep=": "))) +
  geom_point(aes(fill=A.afe.median.mean,size=A.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.A,limits=c(afe.min.A,afe.max.A)) +
  geom_point(data=dat[which(dat$annotation.source=="METABOLIC" & dat$header == "C1 metabolism"),],aes(x="Angelo",y=paste(gene_description,module,sep=": "),color=A.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/A.c1-metabolism.png",width=7.5,height=5.5)

# Hopland
p1<-ggplot(dat[which(dat$annotation.source=="METABOLIC" & dat$header == "C1 metabolism"),],aes(x="Hopland",y=paste(gene_description,module,sep=": "))) +
  geom_point(aes(fill=H.afe.median.mean,size=H.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.H,limits=c(afe.min.H,afe.max.H)) +
  geom_point(data=dat[which(dat$annotation.source=="METABOLIC" & dat$header == "C1 metabolism"),],aes(x="Hopland",y=paste(gene_description,module,sep=": "),color=H.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/H.c1-metabolism.png",width=7.5,height=5.5)

# Sedgwick
p1<-ggplot(dat[which(dat$annotation.source=="METABOLIC" & dat$header == "C1 metabolism"),],aes(x="Sedgwick",y=paste(gene_description,module,sep=": "))) +
  geom_point(aes(fill=S.afe.median.mean,size=S.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.S,limits=c(afe.min.S,afe.max.S)) +
  geom_point(data=dat[which(dat$annotation.source=="METABOLIC" & dat$header == "C1 metabolism"),],aes(x="Sedgwick",y=paste(gene_description,module,sep=": "),color=S.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/S.c1-metabolism.png",width=7.5,height=5.5)

# METABLIC n metabolism
# Angelo
p1<-ggplot(dat[which(dat$annotation.source=="METABOLIC" & dat$header == "Nitrogen cycling"),],aes(x="Angelo",y=paste(gene_description,module,sep=": "))) +
  geom_point(aes(fill=A.afe.median.mean,size=A.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.A,limits=c(afe.min.A,afe.max.A)) +
  geom_point(data=dat[which(dat$annotation.source=="METABOLIC" & dat$header == "Nitrogen cycling"),],aes(x="Angelo",y=paste(gene_description,module,sep=": "),color=A.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/A.n-metabolism.png",width=5.5,height=5.5)

# Hopland
p1<-ggplot(dat[which(dat$annotation.source=="METABOLIC" & dat$header == "Nitrogen cycling"),],aes(x="Hopland",y=paste(gene_description,module,sep=": "))) +
  geom_point(aes(fill=H.afe.median.mean,size=H.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.H,limits=c(afe.min.H,afe.max.H)) +
  geom_point(data=dat[which(dat$annotation.source=="METABOLIC" & dat$header == "Nitrogen cycling"),],aes(x="Hopland",y=paste(gene_description,module,sep=": "),color=H.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/H.n-metabolism.png",width=5.5,height=5.5)

# Sedgwick
p1<-ggplot(dat[which(dat$annotation.source=="METABOLIC" & dat$header == "Nitrogen cycling"),],aes(x="Sedgwick",y=paste(gene_description,module,sep=": "))) +
  geom_point(aes(fill=S.afe.median.mean,size=S.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.S,limits=c(afe.min.S,afe.max.S)) +
  geom_point(data=dat[which(dat$annotation.source=="METABOLIC" & dat$header == "Nitrogen cycling"),],aes(x="Sedgwick",y=paste(gene_description,module,sep=": "),color=S.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/S.n-metabolism.png",width=5.5,height=5.5)

# METABLIC o2 metabolism
# Angelo
p1<-ggplot(dat[which(dat$annotation.source=="METABOLIC" & dat$header == "Oxygen metabolism (Oxidative phosphorylation Complex IV)"),],aes(x="Angelo",y=paste(gene_description,module,sep=": "))) +
  geom_point(aes(fill=A.afe.median.mean,size=A.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.A,limits=c(afe.min.A,afe.max.A)) +
  geom_point(data=dat[which(dat$annotation.source=="METABOLIC" & dat$header == "Oxygen metabolism (Oxidative phosphorylation Complex IV)"),],aes(x="Angelo",y=paste(gene_description,module,sep=": "),color=A.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/A.o2-metabolism.png",width=7,height=5.5)

# Hopland
p1<-ggplot(dat[which(dat$annotation.source=="METABOLIC" & dat$header == "Oxygen metabolism (Oxidative phosphorylation Complex IV)"),],aes(x="Hopland",y=paste(gene_description,module,sep=": "))) +
  geom_point(aes(fill=H.afe.median.mean,size=H.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.H,limits=c(afe.min.H,afe.max.H)) +
  geom_point(data=dat[which(dat$annotation.source=="METABOLIC" & dat$header == "Oxygen metabolism (Oxidative phosphorylation Complex IV)"),],aes(x="Hopland",y=paste(gene_description,module,sep=": "),color=H.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/H.o2-metabolism.png",width=7,height=5.5)

# Sedgwick
p1<-ggplot(dat[which(dat$annotation.source=="METABOLIC" & dat$header == "Oxygen metabolism (Oxidative phosphorylation Complex IV)"),],aes(x="Sedgwick",y=paste(gene_description,module,sep=": "))) +
  geom_point(aes(fill=S.afe.median.mean,size=S.hits),pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=afe.S,limits=c(afe.min.S,afe.max.S)) +
  geom_point(data=dat[which(dat$annotation.source=="METABOLIC" & dat$header == "Oxygen metabolism (Oxidative phosphorylation Complex IV)"),],aes(x="Sedgwick",y=paste(gene_description,module,sep=": "),color=S.sig.ks),shape="*",size=5,show.legend=F,position=position_nudge(0.15)) +
  scale_color_manual(values=c("white","black")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), axis.line = element_blank())

ggsave(p1,file="wsip/results/ks-stats/pathway-bubble-plots/S.o2-metabolism.png",width=7,height=5.5)
