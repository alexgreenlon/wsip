library(vegan)
library(ggplot2)
library(reshape2)
library(plyr)
library(ape)
library(data.table)
library(readxl)
library(dplyr)

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

metacyc.deg<-fread("wsip/results/metacyc/metacyc-pwy-ids-2-degredation-pathways.tsv",header=T,sep='\t',fill=T)
#bin.info<-fread("soil-C-SFA/wsip/binning.d/metawrap.bin_refinement_comp50_cont_25.4.2.20/metawrap-bins.ggkbase-organism_info.drep-set.txt",header=T)
bin.info<-fread("soil-C-SFA/docs/wsip-metawrap-bins.genome-info.drep-filt.csv",sep='\t',header=T)

mycols.tbl<-read.csv("soil-C-SFA/docs/ggkbase_color_scheme.wsip-gtdb-phyla.csv",header=F)
mycols<-as.character(mycols.tbl$V2)
names(mycols)<-as.character(mycols.tbl$V1)


bin.info.pwys<-data.frame(matrix(ncol=(1+length(unique(metacyc.deg$non.redundant))),nrow=0))
colnames(bin.info.pwys)<-c("bin",unique(metacyc.deg$non.redundant))

for (i in  1:length(bin.info$me)){
  metacyc.dat<-fread(paste("wsip/results/metacyc/metawrap-bins/",bin.info$me[i],".report",sep=""),header=F)
  deg.counts<-as.data.frame(merge(metacyc.dat[V8==1],metacyc.deg,by.x="V14",by.y="pathway.code")[,.N,by=list(non.redundant)])
  deg.counts.t<-setNames(data.frame(t(deg.counts[,-1])),deg.counts[,1])
  cbind(data.frame(bin=c(bin.info$me[i])),deg.counts.t)
  bin.info.pwys<-bind_rows(bin.info.pwys,deg.counts.t)
}

bin.info$deg.pwys<-as.numeric()
for (i in 1:length(bin.info$me)){
  metacyc.dat<-fread(paste("wsip/results/metacyc/metawrap-bins/",bin.info$me[i],".report",sep=""),header=F)
  bin.info[i]$deg.pwys<-length(unique(merge(metacyc.dat[V8==1],metacyc.deg,by.x="V14",by.y="pathway.code")$non.redundant))
}

ggplot(bin.info,aes(x=afe.median,y=deg.pwys,color=gtdb.phylum,shape=site)) +
  geom_point() + scale_color_manual(values=mycols) + labs(x="Atom Fraction Excess",y="# of MetaCyc pathways labeled \"Degredation\"") + theme_classic() + theme(,text=element_text(size=20))

bin.paths.dist<-vegdist(bin.info.pwys)
### or
bin.paths.dist<-vegdist(bin.info.pwys)

bin.paths.pcoa<-pcoa(bin.paths.dist)
bin.paths.pcoa.df<-data.frame(bin.paths.pcoa$vectors[,1:2])

bin.paths.pcoa.df$site<-substr(row.names(bin.paths.pcoa.df),1,1)

bin.paths.pcoa.varexp.1<-bin.paths.pcoa$values[1]$Eigenvalues[1]/sum(bin.paths.pcoa$values[1]$Eigenvalues)
bin.paths.pcoa.varexp.2<-bin.paths.pcoa$values[1]$Eigenvalues[2]/sum(bin.paths.pcoa$values[1]$Eigenvalues)

ppp <- ggplot() + coord_fixed() +
  labs(x="Comp1, Axis1", y="Comp2, Axis2") +
  geom_hline(yintercept=0, col="darkgrey") +
  geom_vline(xintercept=0, col="darkgrey") +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", size = 1, linetype = "solid"))

## pcoa's
myxlab<-paste("PC1, ",round(100*bin.paths.pcoa.varexp.1,1),"% Variation Explained",sep="")
myylab<-paste("PC2, ",round(100*bin.paths.pcoa.varexp.2,1),"% Variation Explained",sep="")

ppp + geom_point(data=bin.paths.pcoa.df,aes(x=Axis.1,y=Axis.2)) + xlab(myxlab) + ylab(myylab)

## pcoa colored by site
ppp + geom_point(data=bin.paths.pcoa.df,aes(x=Axis.1,y=Axis.2,color=site)) + xlab(myxlab) + ylab(myylab)

### pcoa colored by activity
ppp + geom_point(data=bin.paths.pcoa.df,aes(x=Axis.1,y=Axis.2,color=bin.info$afe.median)) + xlab(myxlab) + ylab(myylab)

### pcoa colored by phylum
ppp + geom_point(data=bin.paths.pcoa.df,aes(x=Axis.1,y=Axis.2,color=bin.info$phylum)) + xlab(myxlab) + ylab(myylab)



cor.test(x=bin.info$afe.median,y=bin.info$deg.pwys,use="complete.obs")

# 	Pearson's product-moment correlation
#
# data:  bin.info$afe.median and bin.info$deg.pwys
# t = 2.1145, df = 405, p-value = 0.03508
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.007367416 0.199672331
# sample estimates:
#       cor
# 0.1044965

bin.info$site<-substr(bin.info$name,1,1)
for (i in 1:length(unique(bin.info$site))){
  the.site<-unique(bin.info$site)[i]
  site.bin.info<-bin.info[site==the.site]
  if(dim(site.bin.info)[1]>2){
    print(the.site)
    print(cor.test(x=site.bin.info$afe.median,y=site.bin.info$deg.pwys,use="complete.obs"))
  }
}

# [1] "H"
#
# 	Pearson's product-moment correlation
#
# data:  site.bin.info$afe.median and site.bin.info$deg.pwys
# t = 0.42997, df = 191, p-value = 0.6677
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1106295  0.1715831
# sample estimates:
#        cor
# 0.03109653
#
# [1] "A"
#
# 	Pearson's product-moment correlation
#
# data:  site.bin.info$afe.median and site.bin.info$deg.pwys
# t = 2.3195, df = 96, p-value = 0.02249
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.0334706 0.4100399
# sample estimates:
#       cor
# 0.2303615
#
# [1] "S"
#
# 	Pearson's product-moment correlation
#
# data:  site.bin.info$afe.median and site.bin.info$deg.pwys
# t = 1.2546, df = 114, p-value = 0.2122
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.0670452  0.2927852
# sample estimates:
#       cor
# 0.1166979


for (i in 1:length(unique(bin.info$phylum))){
  the.phylum<-unique(bin.info$phylum)[i]
  phylum.bin.info<-bin.info[phylum==the.phylum]
  if(dim(phylum.bin.info)[1]>2){
    print(the.phylum)
    print(cor.test(x=phylum.bin.info$afe.median,y=phylum.bin.info$deg.pwys,use="complete.obs"))
  }
}
#
# [1] "Acidobacteria"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$deg.pwys
# t = -0.28883, df = 17, p-value = 0.7762
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.5079671  0.3969270
# sample estimates:
#         cor
# -0.06988102
#
# [1] "Actinobacteria"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$deg.pwys
# t = 1.3982, df = 174, p-value = 0.1638
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.04318657  0.24943774
# sample estimates:
#       cor
# 0.1054069
#
# [1] "Bacteroidetes"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$deg.pwys
# t = -0.64604, df = 7, p-value = 0.5388
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.7786639  0.5067420
# sample estimates:
#      cor
# -0.23721
#
# [1] "Bdellovibrio"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$deg.pwys
# t = -0.18968, df = 1, p-value = 0.8807
# alternative hypothesis: true correlation is not equal to 0
# sample estimates:
#        cor
# -0.1863615
#
# [1] "Chloroflexi"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$deg.pwys
# t = -0.32706, df = 7, p-value = 0.7532
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.7275277  0.5894671
# sample estimates:
#        cor
# -0.1226818
#
# [1] "Elusimicrobia"
# Error in cor.test.default(x = phylum.bin.info$afe.median, y = phylum.bin.info$deg.pwys,  :
#   not enough finite observations
# > dim(phylum.bin.info)
# [1]  1 41
# > dim(bin.info)
# [1] 433  41
# > for (i in 1:length(unique(bin.info$phylum))){
# +   the.phylum<-unique(bin.info$phylum)[i]
# +   phylum.bin.info<-bin.info[phylum==the.phylum]
# +   if(dim(phylum.bin.info)[1]>2){
# +     print(the.phylum)
# +     print(cor.test(x=phylum.bin.info$afe.median,y=phylum.bin.info$deg.pwys,use="complete.obs"))
# +   }
# + }
# [1] "Acidobacteria"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$deg.pwys
# t = -0.28883, df = 17, p-value = 0.7762
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.5079671  0.3969270
# sample estimates:
#         cor
# -0.06988102
#
# [1] "Actinobacteria"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$deg.pwys
# t = 1.3982, df = 174, p-value = 0.1638
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.04318657  0.24943774
# sample estimates:
#       cor
# 0.1054069
#
# [1] "Bacteroidetes"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$deg.pwys
# t = -0.64604, df = 7, p-value = 0.5388
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.7786639  0.5067420
# sample estimates:
#      cor
# -0.23721
#
# [1] "Bdellovibrio"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$deg.pwys
# t = -0.18968, df = 1, p-value = 0.8807
# alternative hypothesis: true correlation is not equal to 0
# sample estimates:
#        cor
# -0.1863615
#
# [1] "Chloroflexi"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$deg.pwys
# t = -0.32706, df = 7, p-value = 0.7532
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.7275277  0.5894671
# sample estimates:
#        cor
# -0.1226818
#
# [1] "Gemmatimonadetes"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$deg.pwys
# t = -1.9762, df = 19, p-value = 0.06283
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.71682792  0.02284247
# sample estimates:
#       cor
# -0.412916
#
# [1] "Proteobacteria"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$deg.pwys
# t = 1.6, df = 93, p-value = 0.113
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.03915611  0.35355831
# sample estimates:
#       cor
# 0.1636786
#
# [1] "RIF-CHLX"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$deg.pwys
# t = -2.5205, df = 2, p-value = 0.1279
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.9972925  0.5498448
# sample estimates:
#        cor
# -0.8721051
#
# [1] "Thaumarchaeota"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$deg.pwys
# t = -0.98326, df = 4, p-value = 0.3811
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.9224623  0.5769461
# sample estimates:
#        cor
# -0.4411944
#
# [1] "Verrucomicrobia"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$deg.pwys
# t = -0.43049, df = 2, p-value = 0.7088
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.9784501  0.9302280
# sample estimates:
#        cor
# -0.2912064
#
# [1] "unknown"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$deg.pwys
# t = -1.5031, df = 59, p-value = 0.1381
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.42339040  0.06281321
# sample estimates:
#        cor
# -0.1920452
#

gemma.bin.info<-bin.info[phylum=="Gemmatimonadetes"]
ggplot(gemma.bin.info,aes(x=afe.median,y=deg.pwys,color=phylum))+geom_point()

bin.info$carbo.pwys<-as.numeric()
for (i in 1:length(bin.info$name)){
  metacyc.dat<-fread(paste("wsip/results/metacyc/metawrap-bins/",bin.info$name[i],".report",sep=""),header=F)
  metacyc.dat<-merge(metacyc.dat,metacyc.deg,by.x="V14",by.y="pathway.code",all.x=T)
  bin.info[i]$carbo.pwys<-length(unique(metacyc.dat[V8==1 & compound.class == "Carbohydrate Degradation"]$non.redundant))
}

ggplot(bin.info,aes(x=afe.median,y=carbo.pwys,color=gtdb.phylum,shape=site))+geom_point()

cor.test(x=bin.info$afe.median,y=bin.info$carbo.pwys,use="complete.obs")

# Pearson's product-moment correlation
#
# data:  bin.info$afe.median and bin.info$carbo.pwys
# t = 1.4261, df = 405, p-value = 0.1546
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.02670014  0.16674551
# sample estimates:
#      cor
# 0.07068728
#

for (i in 1:length(unique(bin.info$phylum))){
  the.phylum<-unique(bin.info$phylum)[i]
  phylum.bin.info<-bin.info[phylum==the.phylum]
  if(dim(phylum.bin.info)[1]>2){
    print(the.phylum)
    print(cor.test(x=phylum.bin.info$afe.median,y=phylum.bin.info$carbo.pwys,use="complete.obs"))
  }
}

# [1] "Acidobacteria"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$carbo.pwys
# t = -0.942, df = 17, p-value = 0.3594
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.6147476  0.2575330
# sample estimates:
#        cor
# -0.2227298
#
# [1] "Actinobacteria"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$carbo.pwys
# t = 0.88041, df = 174, p-value = 0.3799
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.08213339  0.21242328
# sample estimates:
#        cor
# 0.06659576
#
# [1] "Bacteroidetes"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$carbo.pwys
# t = -0.51775, df = 7, p-value = 0.6206
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.7593229  0.5410859
# sample estimates:
#        cor
# -0.1920475
#
# [1] "Bdellovibrio"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$carbo.pwys
# t = NA, df = 1, p-value = NA
# alternative hypothesis: true correlation is not equal to 0
# sample estimates:
# cor
#  NA
#
# [1] "Chloroflexi"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$carbo.pwys
# t = -0.7725, df = 7, p-value = 0.4651
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.7961949  0.4716366
# sample estimates:
#        cor
# -0.2802735
#
# [1] "Gemmatimonadetes"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$carbo.pwys
# t = -1.0898, df = 19, p-value = 0.2894
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.6103373  0.2112479
# sample estimates:
#        cor
# -0.2425585
#
# [1] "Proteobacteria"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$carbo.pwys
# t = 4.8522, df = 93, p-value = 4.897e-06
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.2726183 0.5969344
# sample estimates:
#       cor
# 0.4494657
#
# [1] "RIF-CHLX"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$carbo.pwys
# t = -0.72101, df = 2, p-value = 0.5458
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.9852155  0.8995727
# sample estimates:
#        cor
# -0.4542078
#
# [1] "Thaumarchaeota"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$carbo.pwys
# t = -2.4833, df = 4, p-value = 0.06797
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.97446309  0.08897213
# sample estimates:
#        cor
# -0.7788253
#
# [1] "Verrucomicrobia"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$carbo.pwys
# t = 0.98831, df = 2, p-value = 0.4272
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.8638313  0.9892794
# sample estimates:
#       cor
# 0.5728244
#
# [1] "unknown"
#
# 	Pearson's product-moment correlation
#
# data:  phylum.bin.info$afe.median and phylum.bin.info$carbo.pwys
# t = -1.0856, df = 59, p-value = 0.2821
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.3784219  0.1159711
# sample estimates:
#        cor
# -0.1399363
#
# Warning message:
# In cor(x, y) : the standard deviation is zero

gemma.bin.info<-bin.info[phylum=="Proteobacteria"]
ggplot(gemma.bin.info,aes(x=afe.median,y=deg.pwys,color=phylum))+geom_point()
