require(vegan)
require(tibble)
require(stringr)
#require(data.table)
require(reshape2)

normcov.trans<-read.csv("soil-C-SFA/docs/wsip-metawrap-bins.genome-info.drep-filt.csv",sep='\t',header=T)
normcov.trans.melt<-melt(normcov.trans, id.vars = c("sample","me","gtdb.phylum"),
                measure.vars = c("rel.abund"))
# the relative abundance numbers above only show the relaive abundance of a bin in the sample
# it was assembled from. to compare sites. We also mapped reads from each sample to all bins (at 96.5% ani)
# but it looks like the samples from each site still separate almost entirely, so seems better to
# calculate community dissimilarities at a higher level, e.g. genus

# a lot of the bins weren't assigned to a genus by gtdb so instead using  an aggregate of all gtdb taxonomy levels

normcov.trans$gtdb<-paste(normcov.trans$gtdb.domain,";",normcov.trans$gtdb.phylum,";",normcov.trans$gtdb.class,";",normcov.trans$gtdb.order,";",normcov.trans$gtdb.family,";",normcov.trans$gtdb.genus,";",normcov.trans$gtdb.species)

# ordination analysis based on relative abundance
rel.abund.long<-aggregate(normcov.trans$rel.abund, by=list(sample=normcov.trans$sample,gtdb=normcov.trans$gtdb),FUN=sum)
rel.abund.wide<-dcast(rel.abund.long,sample ~ gtdb)
rel.abund.wide[is.na(rel.abund.wide)]=0

rownames(rel.abund.wide)<-rel.abund.wide[,1]
rel.abund.wide[,1]<-NULL

rel.abund.dist<-vegdist(rel.abund.wide,method="bray")

site.dat<-read.csv("/Users/alexgreenlon/soil-C-SFA/docs/Greenlon-et-al.2021.supplemental-Supp.table.1.csv",header=T)
site.dat.sub<-site.dat[,c(1,7:11)]
rownames(site.dat.sub)<-site.dat.sub[,1]
site.dat.sub[,1]<-NULL


cap0<-capscale(rel.abund.dist~1,data=site.dat.sub)
cap1<-capscale(rel.abund.dist~.,data=site.dat.sub)
step.res<-ordistep(cap0,scope=formula(cap1),direction="both",perm.max=9999)

step.res<-ordistep(cap0,scope=formula(cap1),direction="both",perm.max=9999)

# Start: rel.abund.dist ~ 1
#
#                                     Df    AIC      F Pr(>F)
# + Avg.MAP1                           1 29.805 2.7163  0.045 *
# + pH                                 1 30.187 2.3234  0.055 .
# + theta.pwp                          1 31.105 1.4129  0.170
# + Moisture....during.incubation.     1 31.287 1.2374  0.240
# + Field.moisture.....at.collection.  1 31.378 1.1509  0.275
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Step: rel.abund.dist ~ Avg.MAP1
#
#            Df    AIC      F Pr(>F)
# - Avg.MAP1  1 30.628 2.7163  0.025 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#                                     Df   AIC     F Pr(>F)
# + theta.pwp                          1 30.33 1.281  0.235
# + Moisture....during.incubation.     1 30.33 1.281  0.245
# + Field.moisture.....at.collection.  1 30.33 1.281  0.255
# + pH                                 1 30.33 1.281  0.280

r2.step.res<-ordiR2step(cap0,scope=formula(cap1),direction="both",perm.max=9999)

# ordination analysis based on activity
afe.long<-aggregate(normcov.trans$afe.median, by=list(sample=normcov.trans$sample,gtdb=normcov.trans$gtdb),FUN=mean)
afe.wide<-dcast(afe.long,sample ~ gtdb)
rownames(afe.wide)<-afe.wide[,1]
afe.wide[,1]<-NULL
afe.wide<-afe.wide+0.1642
afe.wide[is.na(afe.wide)]=0

afe.dist<-vegdist(afe.wide,method="bray")

afe.cap0<-capscale(rel.abund.dist~1,data=site.dat.sub)
afe.cap1<-capscale(rel.abund.dist~.,data=site.dat.sub)
afe.step.res<-ordistep(cap0,scope=formula(cap1),direction="forward",perm.max=9999)

#
# > afe.step.res<-ordistep(cap0,scope=formula(cap1),direction="forward",perm.max=9999)
#
# Start: rel.abund.dist ~ 1
#
#                                     Df    AIC      F Pr(>F)
# + Avg.MAP1                           1 29.805 2.7163  0.055 .
# + pH                                 1 30.187 2.3234  0.065 .
# + theta.pwp                          1 31.105 1.4129  0.195
# + Moisture....during.incubation.     1 31.287 1.2374  0.230
# + Field.moisture.....at.collection.  1 31.378 1.1509  0.285
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
