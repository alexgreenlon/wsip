library(data.table)
library(ggplot2)

afe<-read.csv("/Users/alexgreenlon/soil-C-SFA/wsip/results/AFE/metawrap-qsip-comparisons.csv",header=T,sep=",")
mydat<-read.csv("soil-C-SFA/docs/wsip-metawrap-bins.genome-info.10.7.20.csv",header=T,sep='\t')
afe.meta<-merge(afe,mydat,by.x="bin",by.y="me")

# #phage<-fread("wsip/results/phage/SFA_bulk_metagenomes_virus_host_predictions.tsv")
# phage[,c("domain","phylum","class","order","family","genus","species"):=tstrsplit(GTDB.taxonomy,";",fixed=T)]
# phage$phage.site<-substr(phage$Contig,1,1)
# phage$host.site<-substr(phage$bin.short.name,1,1)
#
# phage.filt<-phage[(is.na(phage$Prophage.blast.score)),]
# phage.filt<-phage.filt[phage.site==host.site,]
# phage.filt<-phage.filt[phage.filt[,.I[Final.score == max(Final.score)],by=Viral.population]$V1]
#
# phage.filt.high<-phage.filt[Final.score>1.5,]

# mycols.tbl<-read.csv("soil-C-SFA/docs/ggkbase_color_scheme.uppercase.csv",header=F)
# mycols<-as.character(mycols.tbl$V2)
# names(mycols)<-as.character(mycols.tbl$V1)

mycols.tbl<-read.csv("soil-C-SFA/docs/ggkbase_color_scheme.wsip-gtdb-phyla.csv",header=F)
mycols<-as.character(mycols.tbl$V2)
names(mycols)<-as.character(mycols.tbl$V1)

#fig.4<-ggplot(afe.16s,aes(x=bin.afe.median,y=X16S.afe.median,shape=site,color=gtdb.phylum))+geom_point() + scale_color_manual(values=mycols) +

afe.meta$AFE.16O<-as.numeric(afe.meta$AFE.16O)
afe.meta$AFE.18O_only.obs<-as.numeric(afe.meta$AFE.18O_only.obs)

coef(lm(AFE.18O_only.obs~AFE.16O,data=afe.meta))
cor.test(afe.meta$AFE.16O,afe.meta$AFE.18O_only.obs)

fig.s8<-ggplot(afe.meta[which(afe.meta$drep.set==1),],aes(x=AFE.16O,y=AFE.18O_only.obs,shape=site,color=gtdb.phylum))+geom_point() + scale_color_manual(values=mycols) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position = "none") +
  xlab("AFE calculated as per Hungate et al., 2015") +  ylab("AFE calculated using only 18O samples")

afe.meta$afe_18O_only_error<-afe.meta$AFE.18O_only.obs-afe.meta$AFE.16O

fig.s8b<-ggplot(afe.meta[which(afe.meta$drep.set==1),],aes(x=afe_18O_only_error,y=gc_avg,shape=site,color=gtdb.phylum))+geom_point() + scale_color_manual(values=mycols) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position = "none") +
  xlab("(18O-only AFE) - (AFE calculated as per Hungate 2015)") +  ylab("GC") + geom_vline(xintercept=0)

# ggplot(phage.filt,aes(x=host.AFE,y=phage.AFE,shape=host.site,col=host.site,size=host.RA))+geom_point()
# ggplot(phage.filt,aes(x=host.AFE,y=phage.AFE,shape=host.site,col=host.RA,size=host.RA))+geom_point()
# ggplot(phage.filt,aes(x=host.AFE,y=phage.AFE,shape=host.site,col=ggkbase.phylum,size=host.RA))+geom_point() + scale_color_manual(values=mycols)
# ggplot(phage.filt[Final.score>1.5,],aes(x=host.AFE,y=phage.AFE,shape=host.site,col=host.site,size=host.RA))+geom_point()
# ggplot(phage.filt[Final.score>1.5,],aes(x=host.AFE,y=phage.AFE,shape=host.site,col=ggkbase.phylum,size=host.RA))+geom_point() + scale_color_manual(values=mycols)
# ggplot(phage.filt[!(is.na(CRISPR.score)),],aes(x=host.AFE,y=phage.AFE,shape=host.site,col=host.site,size=host.RA))+geom_point()
#
#
# ggsave("/Users/alexgreenlon/soil-C-SFA/wsip/results/phage/fig-4.mag-vs-phage-afe.phylum-colors.crispr-matches-only.png",fig.4)
