library(data.table)
library(ggplot2)

ssu.vs.mag<-read.csv("/Users/alexgreenlon/soil-C-SFA/wsip/results/16S-comparison/all.16s-all-contigs.bs.txt",sep="\t",header=T)

phage<-fread("wsip/results/phage/SFA_bulk_metagenomes_virus_host_predictions.tsv")
phage[,c("domain","phylum","class","order","family","genus","species"):=tstrsplit(GTDB.taxonomy,";",fixed=T)]
phage$phage.site<-substr(phage$Contig,1,1)
phage$host.site<-substr(phage$bin.short.name,1,1)

phage.filt<-phage[(is.na(phage$Prophage.blast.score)),]
phage.filt<-phage.filt[phage.site==host.site,]
phage.filt<-phage.filt[phage.filt[,.I[Final.score == max(Final.score)],by=Viral.population]$V1]

phage.filt.high<-phage.filt[Final.score>1.5,]

mycols.tbl<-read.csv("soil-C-SFA/docs/ggkbase_color_scheme.uppercase.csv",header=F)
mycols<-as.character(mycols.tbl$V2)
names(mycols)<-as.character(mycols.tbl$V1)

ggplot(phage.filt,aes(x=host.AFE,y=phage.AFE,shape=host.site,col=host.site,size=host.RA))+geom_point()
ggplot(phage.filt,aes(x=host.AFE,y=phage.AFE,shape=host.site,col=host.RA,size=host.RA))+geom_point()
ggplot(phage.filt,aes(x=host.AFE,y=phage.AFE,shape=host.site,col=ggkbase.phylum,size=host.RA))+geom_point() + scale_color_manual(values=mycols)
ggplot(phage.filt[Final.score>1.5,],aes(x=host.AFE,y=phage.AFE,shape=host.site,col=host.site,size=host.RA))+geom_point()
ggplot(phage.filt[Final.score>1.5,],aes(x=host.AFE,y=phage.AFE,shape=host.site,col=ggkbase.phylum,size=host.RA))+geom_point() + scale_color_manual(values=mycols)
ggplot(phage.filt[!(is.na(CRISPR.score)),],aes(x=host.AFE,y=phage.AFE,shape=host.site,col=host.site,size=host.RA))+geom_point()
fig.4<-ggplot(phage.filt[!(is.na(CRISPR.score)),],aes(x=host.AFE,y=phage.AFE,shape=host.site,col=ggkbase.phylum,size=host.RA))+geom_point() + scale_color_manual(values=mycols) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Host AFE") +  ylab("Phage AFE")

ggsave("/Users/alexgreenlon/soil-C-SFA/wsip/results/phage/fig-4.mag-vs-phage-afe.phylum-colors.crispr-matches-only.png",fig.4)

cor.test(phage.filt[!(is.na(CRISPR.score)),]$host.AFE,phage.filt[!(is.na(CRISPR.score)),]$phage.AFE)
