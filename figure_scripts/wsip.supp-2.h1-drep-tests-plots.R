library(ggplot2)
library("wesanderson")

## plots for the distribution of n50 values for each genome (from each co-assembly) in clusters that had several bins.
# load data
h16.top<-read.csv("soil-C-SFA/wsip/binning.d/H1-16-all-co-assemblies/h1-16.high-binned-drep-genomes.info.csv",header=T)
h18.top<-read.csv("soil-C-SFA/wsip/binning.d/H1-18-all-co-assemblies/h1-18.high-binned-drep-genomes.info.csv",header=T)

# H1-16 plots
p<-ggplot(h16.top,aes(x=Density,y=N50))+geom_bar(stat="identity",colour="black",fill="white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(limits=c("1.6900-1.7720 (all)","1.7350-1.7720","1.7300-1.7468","1.7250-1.7399","1.7200-1.7349","1.7150-1.7299","1.7100-1.7249","1.6900-1.7199")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "blue"))

p + facet_wrap(~secondary_cluster,nrow=1,scales="free_y") +
  xlab("Density of fractions included in assembly") +
  ggtitle("N50's of genomes from each Hopland (rep 1, 16O) co-assembly for most-binned genome ANI clusters")

ggsave("soil-C-SFA/wsip/binning.d/H1-16-all-co-assemblies/h1-16-density-assemblies-vs-n50.png",width=10.2,height=2.88,units="in")

# H1-18 plots
p18<-ggplot(h18.top,aes(x=Density,y=N50))+geom_bar(stat="identity",colour="black",fill="white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(limits=c("1.6900-1.7720 (all)","1.7350-1.7720","1.7300-1.7468","1.7250-1.7399","1.7200-1.7349","1.7150-1.7299","1.7100-1.7249","1.6900-1.7199")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "red"))

p18 + facet_wrap(~secondary_cluster,nrow=1,scales="free_y") +
  xlab("Density of fractions included in assembly") +
  ggtitle("N50's of genomes from each Hopland (rep 1, 18O) co-assembly for most-binned genome ANI clusters")

ggsave("soil-C-SFA/wsip/binning.d/H1-18-all-co-assemblies/h1-18-density-assemblies-vs-n50.png",width=10.2,height=2.88,units="in")

## plots for the number of genomes in each cluster, and which co-assembly the best was from.
# load data
h16.num<-read.csv("soil-C-SFA/wsip/binning.d/H1-16-all-co-assemblies/h1-16.num-genomes-in-clusters.csv",header=T)
h18.num<-read.csv("soil-C-SFA/wsip/binning.d/H1-18-all-co-assemblies/h1-18.num-genomes-in-clusters.csv",header=T)

# make the color scheme (cherry picked from the Moonrise2 and Moonrise3 colors from the wesanders package)
dencolors.wes<-c("1.6900-1.7199"="#C27D38","1.7100-1.7249"="#29211F","1.7150-1.7299"="#85D4E3","1.7200-1.7349"="#F4B5BD","1.7250-1.7399"="#9C964A","1.7300-1.7468"="#CDC08C","1.7350-1.7720"="#FAD77B","1.6900-1.7720 (all)"="#798E87")

ggplot(h16.num,aes(x=reorder(cluster,-genomes.in.cluster),y=genomes.in.cluster,fill=density)) +
  geom_bar(stat="identity")+scale_fill_manual("Density fractions in co-assembly\nof best genome in cluster",values=dencolors.wes,drop=FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "blue")) +
  theme(axis.text.x = element_blank(),axis.ticks.x=element_blank()) +
  xlab("Genome bin ANI cluster") + ylab("# of genome bins in cluster") +
  ggtitle("Number of genome bins from different co-assemblies in ANI clusters from Hopland rep 1, 16O")

ggsave("soil-C-SFA/wsip/binning.d/H1-16-all-co-assemblies/h1-16-num-genomes-vs-cluster.png",width=10.2,height=2.88,units="in")

# note here that i mispelled density in this csv, so needed to change the fill value.
ggplot(h18.num,aes(x=reorder(cluster,-genomes.in.cluster),y=genomes.in.cluster,fill=desnity)) +
  geom_bar(stat="identity") + scale_fill_manual("Density fractions in co-assembly\nof best genome in cluster",values=dencolors.wes,drop=FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "blue")) +
  theme(axis.text.x = element_blank(),axis.ticks.x=element_blank()) +
  xlab("Genome bin ANI cluster") + ylab("# of genome bins in cluster") +
  ggtitle("Number of genome bins from different co-assemblies in ANI clusters from Hopland rep 1, 16O")

ggsave("soil-C-SFA/wsip/binning.d/H1-18-all-co-assemblies/h1-18-num-genomes-vs-cluster.png",width=10.2,height=2.88,units="in")
