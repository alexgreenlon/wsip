library(reshape2)
library(vegan)
library(ape)
library(ggplot2)
library(plyr)
library(data.table)
library(gggenes)

variability_table <- function(cca){
        chi <- c(cca$tot.chi,
                       cca$CCA$tot.chi, cca$CA$tot.chi)
        variability_table <- cbind(chi, chi/chi[1])
        colnames(variability_table) <- c("inertia", "proportion")
        rownames(variability_table) <- c("total", "constrained", "unconstrained")
        return(variability_table)
}
cca_ci <- function(cca, permutations=5000){
        var_tbl <- variability_table(cca)
        p <- permutest(cca, permutations=permutations)
        ci <- quantile(p$F.perm, c(.05,.95),names=F)*p$chi[1]/var_tbl["total", "inertia"]
        return(ci)
}
cap_var_props <- function(cca){
        eig_tot <- sum(cca$CCA$eig)
        var_propdf <- cca$CCA$eig/eig_tot
        return(var_propdf)
}

eps.clusters<-read.table("wsip/results/eps/synth-eps.antismash-clusters/orf2subfamily2family.tsv",header=T,sep='\t')
temp<-do.call(rbind,strsplit(as.character(eps.clusters$orf),split=as.character('|'),fixed=T))
eps.clusters$region<-paste(temp[,1],temp[,2],temp[,3],sep='|')

eps.clusters.info<-read.csv("wsip/results/eps/all-metawrap-bins.synth-eps.antismash.regions.csv",header=T)
# need to figure out what the duplicates are...
#eps.clusters.info<-eps.clusters.info[!duplicated(eps.clusters.info$region),]
metawrap.bin.names<-read.table("wsip/binning.d/metawrap.bin_refinement_comp50_cont_25.4.2.20/metawrap-bins.all-names.txt",header=F,sep='\t')
eps.clusters.info<-merge(eps.clusters.info,data.frame(oldbin=metawrap.bin.names$V2,newbin=paste(metawrap.bin.names$V3,paste("bin.",metawrap.bin.names$V4,sep=""),sep="|")),by.x="bin",by.y="oldbin",keep.x=T)
temp<-do.call(rbind,strsplit(gsub(".gbk","",eps.clusters.info$region),split=as.character('.'),fixed=T))
eps.clusters.info$newregion<-paste(eps.clusters.info$newbin,paste(temp[,3],".",temp[,6],sep="."),sep="|")
eps.clusters.info<-eps.clusters.info[!duplicated(eps.clusters.info$newregion),]
rownames(eps.clusters.info)<-eps.clusters.info$newregion

eps.clusters.wide<-dcast(data=eps.clusters,region~family,value.var="orf")
rownames(eps.clusters.wide)<-eps.clusters.wide$region
eps.clusters.wide<-eps.clusters.wide[,!(names(eps.clusters.wide) %in% c("region"))]
eps.dist<-vegdist(eps.clusters.wide,method="jaccard")

eps.pcoa<-pcoa(eps.dist)

eps.pcoa.df<-data.frame(eps.pcoa$vectors[,1:2])
eps.pcoa.varexp.1<-eps.pcoa$values[1]$Eigenvalues[1]/sum(eps.pcoa$values[1]$Eigenvalues)
eps.pcoa.varexp.2<-eps.pcoa$values[1]$Eigenvalues[2]/sum(eps.pcoa$values[1]$Eigenvalues)

ppp <- ggplot() + coord_fixed() +
  labs(x="Comp1, Axis1", y="Comp2, Axis2") +
  geom_hline(yintercept=0, col="darkgrey") +
  geom_vline(xintercept=0, col="darkgrey") +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", size = 1, linetype = "solid"))

## pcoa's
myxlab<-paste("PC1, ",round(100*eps.pcoa.varexp.1,1),"% Variation Explained",sep="")
myylab<-paste("PC2, ",round(100*eps.pcoa.varexp.2,1),"% Variation Explained",sep="")

ppp + geom_point(data=eps.pcoa.df,aes(x=Axis.1,y=Axis.2)) + xlab(myxlab) + ylab(myylab)

## pcoa colored by polysaccharide
ppp + geom_point(data=eps.pcoa.df,aes(x=Axis.1,y=Axis.2,color=eps.clusters.info$poly)) + xlab(myxlab) + ylab(myylab)

# gene synteny plots

eps.clusters.coords<-read.csv("wsip/results/eps/all-metawrap-bins.synth-eps.antismash.gene-coords.csv",header=T)
eps.clusters.coords$orf<-paste(eps.clusters.coords$region,eps.clusters.coords$gene,sep="|")
eps.clusters<-merge(eps.clusters,eps.clusters.coords,by.x="orf",by.y="orf")
eps.clusters<-merge(eps.clusters,eps.clusters.info,by.x="region.x",by.y="newregion",keep.x=T)

alg.num.fams<-aggregate(eps.clusters[which(eps.clusters$poly=="Alg"),]$family,by=list(unique.values = eps.clusters[which(eps.clusters$poly=="Alg"),]$family),FUN=length)
alg.num.fams.gt2<-subset(alg.num.fams[order(-alg.num.fams$x),],x > 2)

mycols<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff','#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000' )
mycols<-as.character(c(mycols[1:length(alg.num.fams.gt2$unique.values)],rep('#808080',(length(alg.num.fams$unique.values)-length(alg.num.fams.gt2$unique.values)))))
names(mycols)<-as.character(alg.num.fams[order(-alg.num.fams$x),]$unique.values)

eps.clusters.sub<-data.frame(poly=eps.clusters$poly,start=eps.clusters$start.x,end=eps.clusters$end.x,region=eps.clusters$region.x,family=eps.clusters$family)

eps.clusters.sub<-data.frame(poly=eps.clusters$poly,start=eps.clusters$start,end=eps.clusters$end,region=eps.clusters$region,family=eps.clusters$family)

# ggplot(eps.clusters[which(eps.clusters$poly=="Alg"),], aes(xmin = start, xmax = end, y = region.x, fill = family)) +
#   geom_gene_arrow() +
#   facet_wrap(~ region.x, scales = "free", ncol = 1) +
#   scale_fill_manual(name="grp",values=mycols)
#

dummies <- make_alignment_dummies(
  alg.clusters,
  aes(xmin = start, xmax = end, y = region, id = family),
  on = "fam183"
)

ggplot(alg.clusters, aes(xmin = start, xmax = end, y = region, fill = family)) +
  geom_gene_arrow() +
  geom_blank(data = dummies) +
  facet_wrap(~ region, scales = "free", ncol = 1) +
  scale_fill_manual(name="grp",values=mycols) + 
  theme_genes()


ggplot(eps.clusters.sub[which(eps.clusters.sub$poly=="Alg"),], aes(xmin = start, xmax = end, y = region, fill = family)) +
  geom_gene_arrow() +
  facet_wrap(~ region, scales = "free", ncol = 1) +
  scale_fill_manual(name="grp",values=mycols) + 
  theme_genes()
