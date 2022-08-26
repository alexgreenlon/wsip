library(ggplot2)

mydat<-read.csv("soil-C-SFA/docs/wsip-metawrap-bins.genome-info.10.7.20.csv",header=T,sep='\t')

mydat<-mydat[order(mydat$site,mydat$gtdb.phylum,-mydat$afe.median),]
mydat$graph.order<-mydat$graph.order<-c(1:dim(mydat)[1])

#mydat$site<-substr(mydat$name,1,1)

mycols.tbl<-read.csv("soil-C-SFA/docs/ggkbase_color_scheme.wsip-gtdb-phyla.csv",header=F)
mycols<-as.character(mycols.tbl$V2)
names(mycols)<-as.character(mycols.tbl$V1)

p<-ggplot(mydat[which(mydat$drep.set==1),],aes(x=factor(graph.order),afe.median,ymin=lower.ci,ymax=upper.ci,color=gtdb.phylum)) + geom_pointrange(size=0.25) + geom_linerange() + geom_hline(yintercept=0,linetype='dashed',alpha=0.5) + labs(x="bin",y="Atom Fraction Excess") + theme_classic() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) + scale_color_manual(values=mycols)

#ggsave(p,file="soil-C-SFA/bin/mg-qsip-scripts/2A.metawrap-bins.ggkbase-organism_info.drep-set.ci-plot.10.7.20.png",width=11,height=6,units="in")
ggsave(p,file="~/2A.metawrap-bins.ggkbase-organism_info.drep-set.ci-plot.10.7.20.png",width=11,height=6,units="in")

pA<-ggplot(mydat[which(mydat$site=="A" & mydat$drep.set == 1),],aes(x=factor(graph.order),afe.median,ymin=lower.ci,ymax=upper.ci,color=gtdb.phylum)) + geom_pointrange(size=0.25) + geom_linerange() + geom_hline(yintercept=0,linetype='dashed',alpha=0.5) + labs(x="bin",y="Atom Fraction Excess") + theme_classic() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none",text=element_text(size=20)) + scale_color_manual(values=mycols) + ylim(c(-0.4,0.55))

#ggsave(pA,file="soil-C-SFA/bin/mg-qsip-scripts/2A.metawrap-bins.ggkbase-organism_info.drep-set.ci-plot.A.10.7.20.same-y.png",width=3.5,height=6,units="in")
ggsave(pA,file="~/2A.metawrap-bins.ggkbase-organism_info.drep-set.ci-plot.A.10.7.20.same-y.png",width=3.5,height=6,units="in")

pH<-ggplot(mydat[which(mydat$site=="H" & mydat$drep.set == 1),],aes(x=factor(graph.order),afe.median,ymin=lower.ci,ymax=upper.ci,color=gtdb.phylum)) + geom_pointrange(size=0.25) + geom_linerange() + geom_hline(yintercept=0,linetype='dashed',alpha=0.5) + labs(x="bin",y="Atom Fraction Excess") + theme_classic() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none",text=element_text(size=20)) + scale_color_manual(values=mycols) + ylim(c(-0.4,0.55))

#ggsave(pH,file="soil-C-SFA/bin/mg-qsip-scripts/2A.metawrap-bins.ggkbase-organism_info.drep-set.ci-plot.H.10.7.20.same-y.png",width=3.5,height=6,units="in")
ggsave(pH,file="~/2A.metawrap-bins.ggkbase-organism_info.drep-set.ci-plot.H.10.7.20.same-y.png",width=3.5,height=6,units="in")

pS<-ggplot(mydat[which(mydat$site=="S" & mydat$drep.set == 1),],aes(x=factor(graph.order),afe.median,ymin=lower.ci,ymax=upper.ci,color=gtdb.phylum)) + geom_pointrange(size=0.25) + geom_linerange() + geom_hline(yintercept=0,linetype='dashed',alpha=0.5) + labs(x="bin",y="Atom Fraction Excess") + theme_classic() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none",text=element_text(size=20)) + scale_color_manual(values=mycols) + ylim(c(-0.4,0.55))

#ggsave(pS,file="soil-C-SFA/bin/mg-qsip-scripts/2A.metawrap-bins.ggkbase-organism_info.drep-set.ci-plot.S.10.7.20.same-y.png",width=3.5,height=6,units="in")
ggsave(pS,file="~/2A.metawrap-bins.ggkbase-organism_info.drep-set.ci-plot.S.10.7.20.same-y.png",width=3.5,height=6,units="in")
