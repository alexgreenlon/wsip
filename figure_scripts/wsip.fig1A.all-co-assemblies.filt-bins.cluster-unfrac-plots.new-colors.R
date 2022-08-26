library(data.table)
library(ggplot2)

normcov.trans<-fread("soil-C-SFA/docs/wsip-metawrap-bins.genome-info.drep-filt.csv",sep='\t',header=T)
normcov.trans.melt<-melt(normcov.trans, id.vars = c("sample","me","gtdb.phylum"),
                measure.vars = c("rel.abund"))

mycols.tbl<-read.csv("soil-C-SFA/docs/ggkbase_color_scheme.wsip-gtdb-phyla.csv",header=F)
mycols<-as.character(mycols.tbl$V2)
names(mycols)<-as.character(mycols.tbl$V1)

g<-ggplot(normcov.trans.melt,aes(x=sample,y=value))+geom_bar(stat="identity",aes(fill=gtdb.phylum))+scale_fill_manual(name="grp",values=mycols) +
  #theme(axis.title.x=element_blank(),panel.background = element_blank(), axis.line = element_line(size=1,color="black"), axis.text.x = element_text(angle = 90)) +
  theme(text=element_text(size=20),panel.background = element_blank(), axis.line = element_line(size=1,color="black"), axis.text.x = element_text(angle = 90)) +
  xlab("Sample") + ylab("Relative Abundance")

ggsave(g,file="soil-C-SFA/manuscripts/mg-qSIP-wSIP/figures/AG-mg-qSIP-main-figs.1A.bin-rel-abund.png",height=7.5,width=10)

aov(afe.median~site,dat=normcov.trans)

normcov.trans$site.full<-normcov.trans$site
normcov.trans$site.full[normcov.trans$site.full=="A"]<-"Angelo"
normcov.trans$site.full[normcov.trans$site.full=="S"]<-"Sedgwick"
normcov.trans$site.full[normcov.trans$site.full=="H"]<-"Hopland"

gb<-ggplot(normcov.trans,aes(x=site.full,y=afe.median)) + geom_boxplot() +
  theme(text=element_text(size=20),panel.background = element_blank(), axis.line = element_line(size=1,color="black"), axis.text.x = element_text(angle = 90)) +
  xlab("Site") + ylab("Atom Fraction Excess")

site.anova<-aov(afe.median ~ site, data = normcov.trans)
