# compare distributions of afe.median for functions predicted by annotation pipeline DRAM
# to null distribution of afe.median values for all bins using KS test.
# plot the distributions afe.median values for bins at all and each site for each function

library('readxl')
library('reshape2')
library('purrr')
library('data.table')
library(ggplot2)
library(gridExtra)
library(grid)
library(spatstat)

ks_weighted <- function(vector_1,vector_2,weights_1,weights_2){
    F_vec_1 <- ewcdf(vector_1, weights = weights_1, normalise=FALSE)
    F_vec_2 <- ewcdf(vector_2, weights = weights_2, normalise=FALSE)
    xw <- c(vector_1,vector_2)
    d <- max(abs(F_vec_1(xw) - F_vec_2(xw)))

    ## P-VALUE with NORMAL SAMPLE
    # n_vector_1 <- length(vector_1)
    # n_vector_2<- length(vector_2)
    # n <- n_vector_1 * n_vector_2/(n_vector_1 + n_vector_2)

    # P-VALUE EFFECTIVE SAMPLE SIZE as suggested by Monahan
    n_vector_1 <- sum(weights_1)^2/sum(weights_1^2)
    n_vector_2 <- sum(weights_2)^2/sum(weights_2^2)
    n <- n_vector_1 * n_vector_2/(n_vector_1 + n_vector_2)


    pkstwo <- function(x, tol = 1e-06) {
                if (is.numeric(x))
                    x <- as.double(x)
                else stop("argument 'x' must be numeric")
                p <- rep(0, length(x))
                p[is.na(x)] <- NA
                IND <- which(!is.na(x) & (x > 0))
                if (length(IND))
                    p[IND] <- .Call(stats:::C_pKS2, p = x[IND], tol)
                p
            }

    pval <- 1 - pkstwo(sqrt(n) * d)

# this part below I took from the second answer here: https://stackoverflow.com/questions/51339985/pkstwo-error-in-r
# the whole original function is from https://stackoverflow.com/questions/40044375/how-to-calculate-the-kolmogorov-smirnov-statistic-between-two-weighted-samples/40059727
# both answers use code from the ks test itself to calculate the p-value, but seemingly under different assumptions of what type of test it is...
# need to look at the code to figure out which to use...
# and then tidy up the function and actually call it in the loop over all the genes, etc.

# okay, after reading through the code for the ks.test function from base r, it looks like the original way of calculating the p-value
# in this function from the stack overflow post is correct.
# basically, the function calculates p values different ways based on criteria about the data and the type of test.
# the calculation below (commented out) is for a two-sided test (which ours is), with a sufficiently large sample size
# (which ours usually has, but not always. len(x)*len(y) > 10,000) (plus a cororallary of this, which is that it's comparing two numeric vectors),
# AND that there's no "ties" in the data, meaning no data points in the two input vectors that are equivalent. Since one of our inputs is
# a subset of the other, there are by definition ties (actually all of the values are ties) so it always goes to the final calculation
# using C_pKS2, which is what was used in this function from stack overflow. Also, when i do the calculation manually with the ks stat value
# i get the p value calculated that way.

#    p.value <- 1-.Call(stats:::C_pSmirnov2x,
#                       STAT=ks.test.stats,
#                       length(dram.dt.no.na[gene_id== thegene & site=="S" & hits > 0,]$afe.median),
#                       length(bin.info[site=="S"]$afe.median))

    out <- c(KS_Stat=d, P_value=pval)
    return(out)
}

file<-"wsip/results/dram/metawrap-drep-bins/genome-summaries/metabolism_summary.xlsx"
sheets <-excel_sheets(file)
dram.dat <- map_df(sheets, ~ read_excel(file, sheet = .x))
colnames(dram.dat)<-colnames(read_excel(file,1,col_names=T))


dram.df<-reshape2::melt(dram.dat,value.name="hits")

dram.dt<-as.data.table(dram.df)
dram.dt[,"name":=tstrsplit(variable,".",fixed=T)[1]]

################ want to bring in other functions as well. maybe the eps analysis...
# probably metabolic annotations on unbinned contigs...
# i do need to consider how to interface this with dram, since many of them are the same hmm even, so it seems a little redundant...
# and also might mess up the stats?

# load DRAM product summaries.
dram.prod<-read.csv("wsip/results/dram/metawrap-drep-bins/genome-summaries/product.tsv",row.names=1,sep='\t')[,33:99]
dram.prod.hits<-sapply(dram.prod,as.logical)*1
rownames(dram.prod.hits)<-rownames(dram.prod)
dram.prod.df<-reshape2::melt(dram.prod.hits,value.name="hits")
dram.prod.dt<-as.data.table(dram.prod.df)
dram.prod.dt[,"name":=tstrsplit(Var1,".",fixed=T)[1]]

r.len<-length(dram.prod.dt$name)
dram.prod.new<-data.frame(gene_id=character(r.len),gene_description=character(r.len),module=character(r.len),
  header=character(r.len),subheader=character(r.len),variable=factor(r.len),hits=numeric(r.len),name=character(r.len),stringsAsFactors=F)

dram.prod.dt$gene_description<-dram.prod.dt$Var2
dram.prod.new$gene_id<-dram.prod.dt$Var2
dram.prod.new$name<-dram.prod.dt$name
dram.prod.new$hits<-dram.prod.dt$hits
dram.prod.new$header<-"DRAM product"
dram.prod.new$module<-NA
dram.prod.new$subheader<-NA
dram.prod.new$variable<-dram.prod.dt$Var2


# load eps information:
eps<-fread("wsip/results/eps/all-metawrap-bins.synth-eps.hmm.txt",header=F)
anti<-fread("wsip/results/antismash/all-samples.regions-loci.csv",header=F)
dat<-merge(x=eps,y=anti,by.x=c("V2","V3","V4"),by.y=c("V3","V4","V5"),all=F)
write.table(dat,"wsip/results/eps/all-metawrap-bins.synth-eps.antismash.hmm.txt",sep=",")

# make some new variables to keep track of each gene, each region, etc (in case antismash used the same name)
dat$poly<-substr(dat$V6,1,3)
dat$orf<-paste(dat$V4,".",dat$V1.x,sep="")
dat$region<-paste(dat$V2.y,".",dat$V1.y,sep="")

# filter out hits with high evalues
hits<-dat[V8<0.00001]

# filter out hits to multiple hmms
hits.nodups<-hits[,.SD[which.min(V8)],by=list(orf)]

# calculate how many types of polysaccharide loci had hits in each antismash region
polys.by.regions<-dat[,list(nr.poly=length(unique(poly))),by=list(region=region)]
polys.by.regions.filtered<-hits.nodups[,list(nr.poly=length(unique(poly))),by=list(region=region)]

# calculate how many types of polysaccharide loci had hits in each bin
polys.by.bin<-hits.nodups[,list(nr.poly=length(unique(poly))),by=list(bin=V2.y)]

# pick regions that have at least one other gene from the putative eps operon

# probably do it for each polysaccharide.
# for a given polysaccharide
#   for the regions that contain a hit to the synthase for that polysaccharide
#     is the number of loci for that saccharide greater than 1

hits.nodups.noparas<-hits[,.SD[which.max(V9)],by=list(region,V6)]

genes<-c("BcsA_msa_trim","PelF_msa_trim","PgaC_msa_trim","WssB_msa_trim","Alg8_msa_trim")

eps.operons<-data.frame(region=character(),poly=character())
for (i in 1:length(genes)){
  gene<-genes[i]
  thepoly<-substr(gene,1,3)
  regions<-hits.nodups.noparas[V6==gene]$region
  eps.operons<-rbind(eps.operons,data.frame(region=hits.nodups.noparas[region %in% regions,.N,by=.(region,poly)][poly==thepoly & N>1]$region,poly=thepoly))
}

# get all the data together
eps.operons$bin<-do.call(rbind,strsplit(as.character(eps.operons$region),".",fixed=T))[,1]
eps.operons$phylum<-do.call(rbind,strsplit(do.call(rbind,strsplit(as.character(eps.operons$region),".",fixed=T))[,2],"_"))[,1]

# change eps columns to match dram.
r.len<-length(eps.operons$bin)
eps.new<-data.frame(gene_id=character(r.len),gene_description=character(r.len),module=character(r.len),
  header=character(r.len),subheader=character(r.len),variable=factor(r.len),hits=numeric(r.len),name=character(r.len),stringsAsFactors=F)

eps.new$gene_description<-eps.operons$poly
eps.new$gene_id<-eps.operons$poly
eps.new$name<-eps.operons$bin
eps.new$hits<-1
eps.new$header<-"Synthase-dependent EPS"
eps.new$module<-NA
eps.new$subheader<-NA
eps.new$variable<-NA

#######
# add metabolic annotations

samps<-as.vector(read.csv("soil-C-SFA/docs/wsip.all-samples.sorted.list.txt")[,1])
thepath<-paste("soil-C-SFA/wsip/results/metabolic/all-fraction-co-assembly-out/",samps[1],".METABOLIC_result.xlsx",sep="")
dat<-read_excel(thepath,.name_repair='universal')
s<-strsplit(dat$Gn001.Hits, split=",")
newdat<-data.frame(func=rep(dat$Gene.abbreviation,sapply(s,length)),hit=unlist(s))
newdat<-merge(newdat,dat[,1:11],by.y="Gene.abbreviation",by.x="func")
newdat$contig<-sub("(.*)_.*", "\\1",newdat$hit)
metab.all<-newdat

for(i in 2:length(samps)){
  thepath<-paste("soil-C-SFA/wsip/results/metabolic/all-fraction-co-assembly-out/",samps[i],".METABOLIC_result.xlsx",sep="")
  dat<-read_excel(thepath,.name_repair='universal')
  s<-strsplit(dat$Gn001.Hits, split=",")
  newdat<-data.frame(func=rep(dat$Gene.abbreviation,sapply(s,length)),hit=unlist(s))
  newdat<-merge(newdat,dat[,1:11],by.y="Gene.abbreviation",by.x="func")
  newdat$contig<-sub("(.*)_.*", "\\1",newdat$hit)
  metab.all<-rbind(metab.all,newdat)
}

bins<-fread("wsip/binning.d/metawrap.bin_refinement_comp50_cont_25.4.2.20/all-metawrap-bins.binner-names.contigs.txt",header=F)
bin.info<-fread("wsip/binning.d/metawrap.bin_refinement_comp50_cont_25.4.2.20/metawrap-bins.ggkbase-organism_info.drep-set.txt",header=T)
bin.info$site<-substr(bin.info$code,1,1)

metab.bins<-merge(metab.all,bins,by.x="contig",by.y="V1")

r.len<-length(metab.bins$V2)
metab.new<-data.frame(gene_id=character(r.len),gene_description=character(r.len),module=character(r.len),
  header=character(r.len),subheader=character(r.len),variable=factor(r.len),hits=numeric(r.len),name=character(r.len),stringsAsFactors=F)

metab.new$gene_description<-metab.bins$func
metab.new$gene_id<-metab.bins$Hmm.file
metab.new$name<-metab.bins$V2
metab.new$hits<-1
metab.new$header<-metab.bins$Category
metab.new$module<-metab.bins$Function
metab.new$subheader<-metab.bins$Gene.name
metab.new$variable<-metab.bins$V2

######## need to bring in metabolic annotations on unbinned contigs also

contigs<-unique(metab.all$contig)
wmd.path<-paste("soil-C-SFA/wsip/results/AFE/",samps[1],".contig-WMDs.csv",sep="")
temp.wmd<-fread(wmd.path)[contigName.short %in% contigs]

gc.path<-paste("soil-C-SFA/wsip/docs/",samps[1],"-all-fractions_scaffold_min1000.gc.txt",sep="")
contig.gc<-read.csv(gc.path,sep='\t',header=T)
contig.gc$X.Name<-sub(" .*", "\\1",contig.gc$X.Name)

temp.temp.wmd<-merge(temp.wmd,as.data.table(contig.gc)[X.Name %in% contigs],by.x="contigName.short",by.y="X.Name")

wmds<-temp.temp.wmd

for(i in 2:length(samps)){
  wmd.path<-paste("soil-C-SFA/wsip/results/AFE/",samps[i],".contig-WMDs.csv",sep="")
  temp.wmd<-fread(wmd.path)[contigName.short %in% contigs]

  gc.path<-paste("soil-C-SFA/wsip/docs/",samps[i],"-all-fractions_scaffold_min1000.gc.txt",sep="")
  contig.gc<-read.csv(gc.path,sep='\t',header=T)
  contig.gc$X.Name<-sub(" .*", "\\1",contig.gc$X.Name)

  temp.temp.wmd<-merge(temp.wmd,as.data.table(contig.gc)[X.Name %in% contigs],by.x="contigName.short",by.y="X.Name")

  wmds<-rbind(wmds,temp.temp.wmd)
}

metab.wmds<-merge(metab.all,wmds,by.x="contig",by.y="contigName.short")
metab.wmds$AFE<-(((((metab.wmds$mean.18-metab.wmds$mean.16)/metab.wmds$mean.16)+1)*(0.496*metab.wmds$GC+307.691)-(0.496*metab.wmds$GC+307.691))/((12.07747 + (0.496*metab.wmds$GC+307.691))-(0.496*metab.wmds$GC+307.691)))*(1-0.002000429)

r.len<-length(metab.wmds$V1)
metab.unbin.new<-data.frame(gene_id=character(r.len),gene_description=character(r.len),module=character(r.len),
  header=character(r.len),subheader=character(r.len),variable=factor(r.len),hits=numeric(r.len),name=character(r.len),afe=numeric(r.len),stringsAsFactors=F)

metab.unbin.new$gene_description<-metab.wmds$func
metab.unbin.new$gene_id<-metab.wmds$Hmm.file
metab.unbin.new$name<-metab.wmds$contig
metab.unbin.new$hits<-1
metab.unbin.new$header<-metab.wmds$Category
metab.unbin.new$module<-metab.wmds$Function
metab.unbin.new$subheader<-metab.wmds$Gene.name
metab.unbin.new$variable<-metab.wmds$contig
metab.unbin.new$afe<-metab.wmds$AFE

# remove dram annotations that are redundant with metabolic annotations

ko.list<-sub(".hmm","",unique(metab.new$gene_id))
dram.dt.nonred<-dram.dt[!(gene_id %in% ko.list)]

# concatenate all the annotations (for bins) into one data table
all.bin.annot<-rbind(dram.dt.nonred,metab.new,dram.prod.new,eps.new)

#bin.info<-fread("soil-C-SFA/wsip/binning.d/metawrap.bin_refinement_comp50_cont_25.4.2.20/metawrap-bins.ggkbase-organism_info.drep-set.txt",header=T)
#bin.info$site<-substr(bin.info$name,1,1)
#bin.info.all<-fread("soil-C-SFA/docs/wsip.all-metawrap-bins.genome-info.fixed-afe.9.22.20.txt",header=T)
bin.info.all<-fread("soil-C-SFA/docs/wsip-metawrap-bins.genome-info.10.6.20.csv",header=T)
bin.info.w.nas<-bin.info.all[drep.set==1,]

#dram.dt<-merge(dram.dt,bin.info.w.nas,by.x="name",by.y="me")


#dram.dt.no.na<-dram.dt[!is.na(dram.dt$afe.median) & (dram.dt$ci.range<=0.2),]
bin.info<-bin.info.w.nas[!is.na(bin.info.w.nas$afe.median) & (bin.info.w.nas$ci.range<=0.2),]
myvec<-c("me", "afe.median", "rel.abund","site")
bin.info.lite<-bin.info[, ..myvec]

A.ks.test<-data.frame(gene=as.character(),A.p.val=as.numeric(),A.afe.median.mean=as.numeric(),A.p.val.w=as.numeric(),A.afe.median.mean.w=as.numeric(),A.hits=as.numeric())
H.ks.test<-data.frame(gene=as.character(),H.p.val=as.numeric(),H.afe.median.mean=as.numeric(),H.p.val.w=as.numeric(),H.afe.median.mean.w=as.numeric(),H.hits=as.numeric())
S.ks.test<-data.frame(gene=as.character(),S.p.val=as.numeric(),S.afe.median.mean=as.numeric(),S.p.val.w=as.numeric(),S.afe.median.mean.w=as.numeric(),S.hits=as.numeric())

#for (i in 1:length(unique(dram.dt$gene_id))){
#  thegene<-unique(dram.dt$gene_id)[i]
  #if (length(dram.dt.no.na[gene_id== thegene & site=="A" & hits > 0,]$afe.median)>=3){a.ks.test<-rbind(a.ks.test,data.frame(gene=c(thegene),A.p.val=c(ks.test(dram.dt.no.na[gene_id== thegene & site=="A" & hits > 0,]$afe.median,bin.info[site=="A"]$afe.median,na.action=na.omit))$p.value))}
  #if (length(dram.dt.no.na[gene_id== thegene & site=="H" & hits > 0,]$afe.median)>=3){h.ks.test<-rbind(h.ks.test,data.frame(gene=c(thegene),H.p.val=c(ks.test(dram.dt.no.na[gene_id== thegene & site=="H" & hits > 0,]$afe.median,bin.info[site=="H"]$afe.median,na.action=na.omit))$p.value))}

#### okay, so the weighted-ks function above didn't really work (or at least every p value was waaay to low)
# so I'm going to try a ridiculous kluge, which is to just make a new bin.info table where each bin has as
# many rows as it's relative abundance...

## here's the pseudo code:
# make a new data frame with the same columns as bin.info
# for each bin in bin.info
#   n = the bin's relative abundance * 10,000 (rounded to the nearest integer)
#   repeat that bin's row n times in the new data frame.

# bin.info.ra.counts <- data.frame(matrix(ncol = ncol(bin.info), nrow = 0))
# colnames(bin.info.ra.counts) <- colnames(bin.info.lite)
# for (i in 1:length(unique(bin.info.lite$me))){
#   n<-round(10000*bin.info.lite[i,]$rel.abund)
#   for (j in 1:n){
#     bin.info.ra.counts<-rbind(bin.info.ra.counts,bin.info.lite[i,])
#   }
# }


# this will also require re-working the loop below somewhat.
# seems like instead of having a giant table that's every bin with a hit for
# each dram annotation and all of the information about that bin over and over
# i should just have the bin.info tables and at the beginning of each iteration
# subset those tables for the bins that have each function...
# I think the current loop would be unworkable if I make a new table that's weighted by relative abundance.

###################

# A.ks.test<-data.frame(gene=as.character(),A.p.val=as.numeric(),A.afe.median.mean=as.numeric(),A.p.val.w=as.numeric(),A.afe.median.mean.w=as.numeric(),A.hits=as.numeric())
# H.ks.test<-data.frame(gene=as.character(),H.p.val=as.numeric(),H.afe.median.mean=as.numeric(),H.p.val.w=as.numeric(),H.afe.median.mean.w=as.numeric(),H.hits=as.numeric())
# S.ks.test<-data.frame(gene=as.character(),S.p.val=as.numeric(),S.afe.median.mean=as.numeric(),S.p.val.w=as.numeric(),S.afe.median.mean.w=as.numeric(),S.hits=as.numeric())

A.ks.test<-data.frame(gene=as.character(),A.p.val=as.numeric(),A.afe.median.mean=as.numeric(),A.hits=as.numeric())
H.ks.test<-data.frame(gene=as.character(),H.p.val=as.numeric(),H.afe.median.mean=as.numeric(),H.hits=as.numeric())
S.ks.test<-data.frame(gene=as.character(),S.p.val=as.numeric(),S.afe.median.mean=as.numeric(),S.hits=as.numeric())

for (i in 1:length(unique(all.bin.annot$gene_id))){
#for (i in 1:6){
  thegene<-levels(all.bin.annot$gene_id)[i]
  gene.bin.info <- bin.info.lite[me %in% all.bin.annot[gene_id==thegene & hits > 0,]$name]
#  gene.bin.info.ra.counts <- bin.info.ra.counts[me %in% all.bin.annot[gene_id==thegene & hits > 0,]$name]

  gene.bin.info.a<-gene.bin.info[site=="A",]
  #gene.bin.info.ra.counts.a<-gene.bin.info.ra.counts[site=="A",]
  if (length(gene.bin.info.a$afe.median)>=3){
    A.ks.test<-rbind(A.ks.test,data.frame(gene=c(thegene),
      A.p.val=c(ks.test(gene.bin.info.a$afe.median,
        bin.info.lite[site=="A"]$afe.median,na.action=na.omit)$p.value),
      A.afe.median.mean=c(mean(gene.bin.info.a$afe.median)),
#      A.p.val.w=c(ks.test(gene.bin.info.ra.counts.a$afe.median,bin.info.ra.counts[site=="A"]$afe.median)$p.value),
#      A.afe.median.mean.w=c(weighted.mean(gene.bin.info.a$afe.median,
#        gene.bin.info.a$rel.abund)),
      A.hits=c(length(gene.bin.info.a$afe.median))
      ))
  }

  gene.bin.info.a<-gene.bin.info[site=="H",]
#  gene.bin.info.ra.counts.a<-gene.bin.info.ra.counts[site=="H",]
  if (length(gene.bin.info.a$afe.median)>=3){
    H.ks.test<-rbind(H.ks.test,data.frame(gene=c(thegene),
      H.p.val=c(ks.test(gene.bin.info.a$afe.median,
        bin.info.lite[site=="H"]$afe.median,na.action=na.omit)$p.value),
      H.afe.median.mean=c(mean(gene.bin.info.a$afe.median)),
#      H.p.val.w=c(ks.test(gene.bin.info.ra.counts.a$afe.median,bin.info.ra.counts[site=="H"]$afe.median)$p.value),
#      H.afe.median.mean.w=c(weighted.mean(gene.bin.info.a$afe.median,
#        gene.bin.info.a$rel.abund)),
      H.hits=c(length(gene.bin.info.a$afe.median))
      ))
  }

  gene.bin.info.a<-gene.bin.info[site=="S",]
#  gene.bin.info.ra.counts.a<-gene.bin.info.ra.counts[site=="S",]
  if (length(gene.bin.info.a$afe.median)>=3){
    S.ks.test<-rbind(S.ks.test,data.frame(gene=c(thegene),
      S.p.val=c(ks.test(gene.bin.info.a$afe.median,
        bin.info.lite[site=="S"]$afe.median,na.action=na.omit)$p.value),
      S.afe.median.mean=c(mean(gene.bin.info.a$afe.median)),
#      S.p.val.w=c(ks.test(gene.bin.info.ra.counts.a$afe.median,bin.info.ra.counts[site=="S"]$afe.median)$p.value),
#      S.afe.median.mean.w=c(weighted.mean(gene.bin.info.a$afe.median,
#        gene.bin.info.a$rel.abund)),
      S.hits=c(length(gene.bin.info.a$afe.median))
      ))
  }

}

##### ks tests for unbinned contigs with metabolic annotations

A.ks.test$bins<-TRUE
H.ks.test$bins<-TRUE
S.ks.test$bins<-TRUE

A.ks.test.temp<-A.ks.test
H.ks.test.temp<-H.ks.test
S.ks.test.temp<-S.ks.test

rpl6<-fread("wsip/results/SCGs/all-samples.essential-hmm-hits.AFE.csv")[(gene == "Ribosomal_L6" & is.na(AFE)==FALSE),]
rpl6$site<-substr(rpl6$contig,1,1)

metab.unbin.new<-as.data.table(metab.unbin.new)
metab.unbin.new$site<-substr(metab.unbin.new$name,1,1)
metab.unbin.new$gene_id<-paste(metab.unbin.new$gene_id,".all.contigs",sep="")

A.ks.test.unbin<-data.frame(gene=as.character(),A.p.val=as.numeric(),A.afe.median.mean=as.numeric(),A.hits=as.numeric())
H.ks.test.unbin<-data.frame(gene=as.character(),H.p.val=as.numeric(),H.afe.median.mean=as.numeric(),H.hits=as.numeric())
S.ks.test.unbin<-data.frame(gene=as.character(),S.p.val=as.numeric(),S.afe.median.mean=as.numeric(),S.hits=as.numeric())

metab.unbin.new.w.nas<-metab.unbin.new
metab.unbin.new<-metab.unbin.new.w.nas[!is.na(metab.unbin.new.w.nas$afe),]

for (i in 1:length(unique(metab.unbin.new$gene_id))){
#for (i in 1:6){
  thegene<-unique(metab.unbin.new$gene_id)[i]
  gene.bin.info <- metab.unbin.new[gene_id==thegene & hits > 0,]

  gene.bin.info.a<-gene.bin.info[site=="A",]
  if (length(gene.bin.info.a$afe)>=5){
    A.ks.test.unbin<-rbind(A.ks.test.unbin,data.frame(gene=c(thegene),
      A.p.val=c(ks.test(gene.bin.info.a$afe,
        rpl6[site=="A"]$AFE,na.action=na.omit)$p.value),
      A.afe.median.mean=c(mean(gene.bin.info.a$afe)),
      A.hits=c(length(gene.bin.info.a$afe))
      ))
  }

  gene.bin.info.a<-gene.bin.info[site=="H",]
  if (length(gene.bin.info.a$afe)>=5){
    H.ks.test.unbin<-rbind(H.ks.test.unbin,data.frame(gene=c(thegene),
      H.p.val=c(ks.test(gene.bin.info.a$afe,
        rpl6[site=="H"]$AFE,na.action=na.omit)$p.value),
      H.afe.median.mean=c(mean(gene.bin.info.a$afe)),
      H.hits=c(length(gene.bin.info.a$afe))
      ))
  }

  gene.bin.info.a<-gene.bin.info[site=="S",]
  if (length(gene.bin.info.a$afe)>=5){
    S.ks.test.unbin<-rbind(S.ks.test.unbin,data.frame(gene=c(thegene),
      S.p.val=c(ks.test(gene.bin.info.a$afe,
        rpl6[site=="S"]$AFE,na.action=na.omit)$p.value),
      S.afe.median.mean=c(mean(gene.bin.info.a$afe)),
      S.hits=c(length(gene.bin.info.a$afe))
      ))
  }

}

A.ks.test.unbin$bins<-FALSE
H.ks.test.unbin$bins<-FALSE
S.ks.test.unbin$bins<-FALSE

A.ks.test<-rbind(A.ks.test,A.ks.test.unbin)
H.ks.test<-rbind(H.ks.test,H.ks.test.unbin)
S.ks.test<-rbind(S.ks.test,S.ks.test.unbin)

# run FDR multiple-comparison corrections
A.ks.test$A.p.adj<-p.adjust(A.ks.test$A.p.val,method="BH")
H.ks.test$H.p.adj<-p.adjust(H.ks.test$H.p.val,method="BH")
S.ks.test$S.p.adj<-p.adjust(S.ks.test$S.p.val,method="BH")

# A.ks.test$A.p.w.adj<-p.adjust(A.ks.test$A.p.val.w,method="BH")
# H.ks.test$H.p.w.adj<-p.adjust(H.ks.test$H.p.val.w,method="BH")
# S.ks.test$S.p.w.adj<-p.adjust(S.ks.test$S.p.val.w,method="BH")

# consolidate each dataframe for each site into one.
#afe.unbin<-merge(A.ks.test[which(A.ks.test$bins==FALSE),],merge(H.ks.test[which(H.ks.test$bins==FALSE),],S.ks.test[which(S.ks.test$bins==FALSE),],by.x="gene",by.y="gene",all=T),by.x="gene",by.y="gene",all=T)
#afe.unbin<-merge(A.ks.test.unbin,merge(H.ks.test.unbin,S.ks.test.unbin,by.x="gene",by.y="gene",all=T),by.x="gene",by.y="gene",all=T)
#afe.unbin$clust<-NA

# cluster based on average afe for bins with each trait
# afe<-merge(A.ks.test.temp,merge(H.ks.test.temp,S.ks.test.temp,by.x="gene",by.y="gene",all=T),by.x="gene",by.y="gene",all=T)
# afe.no.na<-merge(A.ks.test.temp,merge(H.ks.test.temp,S.ks.test.temp,by.x="gene",by.y="gene"),by.x="gene",by.y="gene")
# afe.df<-data.frame(row.names=afe.no.na$gene,A=afe.no.na$A.afe.median.mean,H=afe.no.na$H.afe.median.mean,S=afe.no.na$S.afe.median.mean)
# afe.dist<-dist(afe.df,method='euclidean')
# afe.clust<-hclust(afe.dist)
# afe.cuts<-cutree(afe.clust,k=10)
# afe.cuts<-data.frame(gene=row.names(afe.df),clust=afe.cuts)
# afe<-merge(afe,afe.cuts,by.x="gene",by.y="gene",all=T)

afe<-merge(A.ks.test,merge(H.ks.test,S.ks.test,by.x="gene",by.y="gene",all=T),by.x="gene",by.y="gene",all=T)
afe.no.na<-merge(A.ks.test,merge(H.ks.test,S.ks.test,by.x="gene",by.y="gene"),by.x="gene",by.y="gene")
afe.df<-data.frame(row.names=afe.no.na$gene,A=afe.no.na$A.afe.median.mean,H=afe.no.na$H.afe.median.mean,S=afe.no.na$S.afe.median.mean)
afe.dist<-dist(afe.df,method='euclidean')
afe.clust<-hclust(afe.dist)
afe.cuts<-cutree(afe.clust,k=10)
afe.cuts<-data.frame(gene=row.names(afe.df),clust=afe.cuts)
afe<-merge(afe,afe.cuts,by.x="gene",by.y="gene",all=T)

# separately cluster for relative-abundance weighted AFE
# afe.no.na<-merge(A.ks.test,merge(H.ks.test,S.ks.test,by.x="gene",by.y="gene"),by.x="gene",by.y="gene")
# afe.df<-data.frame(row.names=afe.no.na$gene,A=afe.no.na$A.afe.median.mean.w,H=afe.no.na$H.afe.median.mean.w,S=afe.no.na$S.afe.median.mean.w)
# afe.dist<-dist(afe.df,method='euclidean')
# afe.clust<-hclust(afe.dist)
# afe.cuts<-cutree(afe.clust,k=10)
# afe.cuts<-data.frame(gene=row.names(afe.df),clust.w=afe.cuts)
# afe<-merge(afe,afe.cuts,by.x="gene",by.y="gene",all=T)

#genes<-as.data.frame(all.bin.annot[,1:5])
#afe<-merge(afe,genes[!duplicated(genes$gene_id),],by.x="gene",by.y="gene_id",all.x=T)

#genes<-as.data.frame(metab.unbin.new[,1:5])
#afe.unbin<-merge(afe.unbin,genes[!duplicated(genes$gene_id),],by.x="gene",by.y="gene_id",all.x=T)

genes<-rbind(as.data.frame(all.bin.annot[,1:5]),as.data.frame(metab.unbin.new[,1:5]))
afe<-merge(afe,genes[!duplicated(genes$gene_id),],by.x="gene",by.y="gene_id",all.x=T)

#afe.all<-merge(A.ks.test,merge(H.ks.test,S.ks.test,by.x="gene",by.y="gene",all=T),by.x="gene",by.y="gene",all=T)
#afe.all<-merge(rbind(afe,afe.unbin),afe.all[,c("gene","A.p.adj","H.p.adj","S.p.adj")],by.x="gene",by.y="gene",all=T)

write.csv(afe, file="wsip/results/ks-stats/combined-annotations.fixed-afe-median.sitewise-ks-tests.fixed.csv")


afe.mean<-mean(bin.info$afe.median)
afe.mean.w<-weighted.mean(bin.info$afe.median,bin.info$rel.abund)
a.afe.min<-min(bin.info[site=="A"]$afe.median)
a.afe.max<-max(bin.info[site=="A"]$afe.median)
a.afe.mean<-mean(bin.info[site=="A"]$afe.median)
a.afe.w.mean<-weighted.mean(bin.info[site=="A"]$afe.median,bin.info[site=="A"]$rel.abund)
h.afe.min<-min(bin.info[site=="H"]$afe.median)
h.afe.max<-max(bin.info[site=="H"]$afe.median)
h.afe.mean<-mean(bin.info[site=="H"]$afe.median)
h.afe.w.mean<-weighted.mean(bin.info[site=="H"]$afe.median,bin.info[site=="H"]$rel.abund)
s.afe.min<-min(bin.info[site=="S"]$afe.median)
s.afe.max<-max(bin.info[site=="S"]$afe.median)
s.afe.mean<-mean(bin.info[site=="S"]$afe.median)
s.afe.w.mean<-weighted.mean(bin.info[site=="S"]$afe.median,bin.info[site=="S"]$rel.abund)

#mycols.tbl<-read.csv("soil-C-SFA/docs/ggkbase_color_scheme.uppercase.csv",header=F)
mycols.tbl<-read.csv("soil-C-SFA/docs/ggkbase_color_scheme.wsip-gtdb-phyla.csv",header=F)
mycols<-as.character(mycols.tbl$V2)
names(mycols)<-as.character(mycols.tbl$V1)

#bin.info.ra.full<-merge(bin.info.ra.counts,bin.info,by.x="me",by.y="me",all.x=T)
#ptest<-ggplot(bin.info.ra.full,aes(x=afe.median.x,fill=phylum))+geom_histogram(binwidth=0.02) +
#scale_fill_manual(name="grp",values=mycols) + xlim(c(-0.2,0.5)) +
#theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))

p1<-ggplot(bin.info,aes(x=afe.median,fill=gtdb.phylum))+geom_histogram(binwidth=0.02) +
scale_fill_manual(name="grp",values=mycols) + ylim(c(0,50)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept = afe.mean, size = 1) +
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))

p2<-ggplot(bin.info[which(bin.info$site=="A"),],aes(x=afe.median,fill=gtdb.phylum)) + geom_histogram(binwidth=0.02) + geom_vline(xintercept=a.afe.mean, size=1, linetype="dashed") +
scale_fill_manual(name="grp",values=mycols) + ylim(c(0,50)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept=afe.mean, size=1) +
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))

p3<-ggplot(bin.info[which(bin.info$site=="H"),],aes(x=afe.median,fill=gtdb.phylum))+geom_histogram(binwidth=0.02) + geom_vline(xintercept=h.afe.mean, size=1, linetype="dashed") +
scale_fill_manual(name="grp",values=mycols) + ylim(c(0,50)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept=afe.mean, size=1) +
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))

p4<-ggplot(bin.info[which(bin.info$site=="S"),],aes(x=afe.median,fill=gtdb.phylum))+geom_histogram(binwidth=0.02) + geom_vline(xintercept=s.afe.mean, size=1, linetype="dashed") +
scale_fill_manual(name="grp",values=mycols) + ylim(c(0,50)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept=afe.mean, size=1) +
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))

#grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,nrow=2)
g<-arrangeGrob(p1,p2,p3,p4,nrow=1)
path<-paste("wsip/results/dram/metawrap-drep-bins/genome-summaries/all-drep-metawrap-bins.fixed-afe.png")
ggsave(g,file=path,width=10,height=3,units="in")

g<-arrangeGrob(p2,p3,p4,nrow=1)
path<-paste("wsip/results/dram/metawrap-drep-bins/genome-summaries/all-drep-metawrap-bins.sites.fixed-afe.png")
ggsave(g,file=path,width=7.5,height=3,units="in")

## AFE distributions weighted by relative abundance
p1<-ggplot(bin.info,aes(x=afe.median, fill=gtdb.phylum))+geom_histogram(aes(weight = rel.abund), binwidth=0.02) +
scale_fill_manual(name="grp",values=mycols) + ylim(c(0,0.4)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept = afe.mean.w, size = 1) +
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))
p2<-ggplot(bin.info[which(bin.info$site=="A"),],aes(x=afe.median,fill=gtdb.phylum))+geom_histogram(aes(weight = rel.abund),binwidth=0.02) + geom_vline(xintercept=a.afe.w.mean, size=1, linetype="dashed") +
scale_fill_manual(name="grp",values=mycols) + ylim(c(0,0.4)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept = afe.mean.w, size = 1) +
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))
p3<-ggplot(bin.info[which(bin.info$site=="H"),],aes(x=afe.median,fill=gtdb.phylum))+geom_histogram(aes(weight = rel.abund),binwidth=0.02) + geom_vline(xintercept=s.afe.w.mean, size=1, linetype="dashed") +
scale_fill_manual(name="grp",values=mycols) + ylim(c(0,0.4)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept = afe.mean.w, size = 1) +
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))
p4<-ggplot(bin.info[which(bin.info$site=="S"),],aes(x=afe.median,fill=gtdb.phylum))+geom_histogram(aes(weight = rel.abund),binwidth=0.02) + geom_vline(xintercept=s.afe.w.mean, size=1, linetype="dashed") +
scale_fill_manual(name="grp",values=mycols) + ylim(c(0,0.4)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept = afe.mean.w, size = 1) +
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))

#grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,nrow=2)
g<-arrangeGrob(p1,p2,p3,p4,nrow=1)
path<-paste("wsip/results/dram/metawrap-drep-bins/genome-summaries/all-drep-metawrap-bins.ra-weighted.fixed-afe.png")
ggsave(g,file=path,width=10,height=3,units="in")

g<-arrangeGrob(p2,p3,p4,nrow=1)
path<-paste("wsip/results/dram/metawrap-drep-bins/genome-summaries/all-drep-metawrap-bins.sites.ra-weighted.fixed-afe.png")
ggsave(g,file=path,width=7.5,height=3,units="in")

for (i in 1:length(unique(all.bin.annot$gene_id))){
  gene<-unique(all.bin.annot$gene_id)[i]
  #gene<-"K02406"
  thebins<-all.bin.annot[which(all.bin.annot$gene_id==gene & all.bin.annot$hits > 0),]$name
  bin.info.gene<-bin.info[which(bin.info$me %in% thebins),]

  afe.mean.g<-mean(bin.info.gene$afe.median)
  afe.mean.w.g<-weighted.mean(bin.info.gene$afe.median,bin.info.gene$rel.abund)
  a.afe.mean.g<-mean(bin.info.gene[site=="A"]$afe.median)
  a.afe.w.mean.g<-weighted.mean(bin.info.gene[site=="A"]$afe.median,bin.info.gene[site=="A"]$rel.abund)
  h.afe.mean.g<-mean(bin.info.gene[site=="H"]$afe.median)
  h.afe.w.mean.g<-weighted.mean(bin.info.gene[site=="H"]$afe.median,bin.info.gene[site=="H"]$rel.abund)
  s.afe.mean.g<-mean(bin.info.gene[site=="S"]$afe.median)
  s.afe.w.mean.g<-weighted.mean(bin.info.gene[site=="S"]$afe.median,bin.info.gene[site=="S"]$rel.abund)

  p1<-ggplot(bin.info.gene,aes(x=afe.median,fill=gtdb.phylum))+geom_histogram(binwidth=0.02) +
  scale_fill_manual(name="grp",values=mycols) + ylim(c(0,50)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept = afe.mean, size = 1) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))
  p2<-ggplot(bin.info.gene[which(bin.info.gene$site=="A"),],aes(x=afe.median,fill=gtdb.phylum)) + geom_histogram(binwidth=0.02) + geom_vline(xintercept=a.afe.mean.g, size=1, linetype="dashed") +
  scale_fill_manual(name="grp",values=mycols) + ylim(c(0,50)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept=a.afe.mean, size=1) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))
  p3<-ggplot(bin.info.gene[which(bin.info.gene$site=="H"),],aes(x=afe.median,fill=gtdb.phylum))+geom_histogram(binwidth=0.02) + geom_vline(xintercept=h.afe.mean.g, size=1, linetype="dashed") +
  scale_fill_manual(name="grp",values=mycols) + ylim(c(0,50)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept=h.afe.mean, size=1) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))
  p4<-ggplot(bin.info.gene[which(bin.info.gene$site=="S"),],aes(x=afe.median,fill=gtdb.phylum))+geom_histogram(binwidth=0.02) + geom_vline(xintercept=s.afe.mean.g, size=1, linetype="dashed") +
  scale_fill_manual(name="grp",values=mycols) + ylim(c(0,50)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept=s.afe.mean, size=1) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))

  g<-arrangeGrob(p2,p3,p4,nrow=1)
  path<-paste("wsip/results/dram/metawrap-drep-bins/genome-summaries/fixed-afe-plots/",eval(gene),".drep-metawrap-bins.png",sep="")
  #path<-paste("wsip/results/dram/metawrap-drep-bins/genome-summaries/fixed-afe-plots/test.drep-metawrap-bins.png",sep="")
  ggsave(g,file=path,width=7.5,height=3,units="in")

  ## AFE distributions weighted by relative abundance
  p1<-ggplot(bin.info.gene,aes(x=afe.median, fill=gtdb.phylum))+geom_histogram(aes(weight = rel.abund), binwidth=0.02) +
  scale_fill_manual(name="grp",values=mycols) + ylim(c(0,0.4)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept = afe.mean.w, size = 1) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))
  p2<-ggplot(bin.info.gene[which(bin.info.gene$site=="A"),],aes(x=afe.median,fill=gtdb.phylum))+geom_histogram(aes(weight = rel.abund),binwidth=0.02) + geom_vline(xintercept=a.afe.w.mean.g, size=1, linetype="dashed") +
  scale_fill_manual(name="grp",values=mycols) + ylim(c(0,0.4)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept = a.afe.w.mean, size = 1) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))
  p3<-ggplot(bin.info.gene[which(bin.info.gene$site=="H"),],aes(x=afe.median,fill=gtdb.phylum))+geom_histogram(aes(weight = rel.abund),binwidth=0.02) + geom_vline(xintercept=s.afe.w.mean.g, size=1, linetype="dashed") +
  scale_fill_manual(name="grp",values=mycols) + ylim(c(0,0.4)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept = h.afe.w.mean, size = 1) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))
  p4<-ggplot(bin.info.gene[which(bin.info.gene$site=="S"),],aes(x=afe.median,fill=gtdb.phylum))+geom_histogram(aes(weight = rel.abund),binwidth=0.02) + geom_vline(xintercept=s.afe.w.mean.g, size=1, linetype="dashed") +
  scale_fill_manual(name="grp",values=mycols) + ylim(c(0,0.4)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept = s.afe.w.mean, size = 1) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))

  g<-arrangeGrob(p2,p3,p4,nrow=1)
  path<-paste("wsip/results/dram/metawrap-drep-bins/genome-summaries/fixed-afe-plots.ra-weighted/",eval(gene),".drep-metawrap-bins.png",sep="")
  #path<-paste("wsip/results/dram/metawrap-drep-bins/genome-summaries/fixed-afe-plots.ra-weighted/.drep-metawrap-bins.png",sep="")
  ggsave(g,file=path,width=7.5,height=3,units="in")
}

for (i in 1:length(unique(eps.new$gene_id))){
  thegene<-unique(eps.new$gene_id)[i]
  #gene<-"K02406"
  thebins<-eps.new[which(eps.new$gene_id==thegene & eps.new$hits > 0),]$name
  bin.info.gene<-bin.info[which(bin.info$me %in% thebins),]

  afe.mean.g<-mean(bin.info.gene$afe.median)
  afe.mean.w.g<-weighted.mean(bin.info.gene$afe.median,bin.info.gene$rel.abund)
  a.afe.mean.g<-mean(bin.info.gene[site=="A"]$afe.median)
  a.afe.w.mean.g<-weighted.mean(bin.info.gene[site=="A"]$afe.median,bin.info.gene[site=="A"]$rel.abund)
  h.afe.mean.g<-mean(bin.info.gene[site=="H"]$afe.median)
  h.afe.w.mean.g<-weighted.mean(bin.info.gene[site=="H"]$afe.median,bin.info.gene[site=="H"]$rel.abund)
  s.afe.mean.g<-mean(bin.info.gene[site=="S"]$afe.median)
  s.afe.w.mean.g<-weighted.mean(bin.info.gene[site=="S"]$afe.median,bin.info.gene[site=="S"]$rel.abund)

  p1<-ggplot(bin.info.gene,aes(x=afe.median,fill=gtdb.phylum))+geom_histogram(binwidth=0.02) +
  scale_fill_manual(name="grp",values=mycols) + ylim(c(0,10)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept = afe.mean, size = 1) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))
  p2<-ggplot(bin.info.gene[which(bin.info.gene$site=="A"),],aes(x=afe.median,fill=gtdb.phylum)) + geom_histogram(binwidth=0.02) + geom_vline(xintercept=a.afe.mean.g, size=1, linetype="dashed") +
  scale_fill_manual(name="grp",values=mycols) + ylim(c(0,10)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept=a.afe.mean, size=1) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))
  p3<-ggplot(bin.info.gene[which(bin.info.gene$site=="H"),],aes(x=afe.median,fill=gtdb.phylum))+geom_histogram(binwidth=0.02) + geom_vline(xintercept=h.afe.mean.g, size=1, linetype="dashed") +
  scale_fill_manual(name="grp",values=mycols) + ylim(c(0,10)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept=h.afe.mean, size=1) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))
  p4<-ggplot(bin.info.gene[which(bin.info.gene$site=="S"),],aes(x=afe.median,fill=gtdb.phylum))+geom_histogram(binwidth=0.02) + geom_vline(xintercept=s.afe.mean.g, size=1, linetype="dashed") +
  scale_fill_manual(name="grp",values=mycols) + ylim(c(0,10)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept=s.afe.mean, size=1) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))

  g<-arrangeGrob(p2,p3,p4,nrow=1)
  path<-paste("wsip/results/dram/metawrap-drep-bins/genome-summaries/fixed-afe-plots/",eval(thegene),".drep-metawrap-bins.ymax10.png",sep="")
  #path<-paste("wsip/results/dram/metawrap-drep-bins/genome-summaries/fixed-afe-plots/test.drep-metawrap-bins.png",sep="")
  ggsave(g,file=path,width=7.5,height=3,units="in")

  ## AFE distributions weighted by relative abundance
  p1<-ggplot(bin.info.gene,aes(x=afe.median, fill=gtdb.phylum))+geom_histogram(aes(weight = rel.abund), binwidth=0.02) +
  scale_fill_manual(name="grp",values=mycols) + ylim(c(0,0.1)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept = afe.mean.w, size = 1) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))
  p2<-ggplot(bin.info.gene[which(bin.info.gene$site=="A"),],aes(x=afe.median,fill=gtdb.phylum))+geom_histogram(aes(weight = rel.abund),binwidth=0.02) + geom_vline(xintercept=a.afe.w.mean.g, size=1, linetype="dashed") +
  scale_fill_manual(name="grp",values=mycols) + ylim(c(0,0.1)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept = a.afe.w.mean, size = 1) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))
  p3<-ggplot(bin.info.gene[which(bin.info.gene$site=="H"),],aes(x=afe.median,fill=gtdb.phylum))+geom_histogram(aes(weight = rel.abund),binwidth=0.02) + geom_vline(xintercept=s.afe.w.mean.g, size=1, linetype="dashed") +
  scale_fill_manual(name="grp",values=mycols) + ylim(c(0,0.1)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept = h.afe.w.mean, size = 1) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))
  p4<-ggplot(bin.info.gene[which(bin.info.gene$site=="S"),],aes(x=afe.median,fill=gtdb.phylum))+geom_histogram(aes(weight = rel.abund),binwidth=0.02) + geom_vline(xintercept=s.afe.w.mean.g, size=1, linetype="dashed") +
  scale_fill_manual(name="grp",values=mycols) + ylim(c(0,0.1)) + xlim(c(-0.2,0.5)) + geom_vline(xintercept = s.afe.w.mean, size = 1) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank(), legend.position = "none", axis.line = element_line(colour = "black"))

  g<-arrangeGrob(p2,p3,p4,nrow=1)
  path<-paste("wsip/results/dram/metawrap-drep-bins/genome-summaries/fixed-afe-plots.ra-weighted/",eval(thegene),".drep-metawrap-bins.ymax0.1.png",sep="")
  #path<-paste("wsip/results/dram/metawrap-drep-bins/genome-summaries/fixed-afe-plots.ra-weighted/.drep-metawrap-bins.png",sep="")
  ggsave(g,file=path,width=7.5,height=3,units="in")
}

mean(bin.info[phylum=="Actinobacteria"]$afe.median)
#gene<-unique(all.bin.annot$gene_id)[i]
gene<-"K03520.hmm"
thebins<-metab.new[which(metab.new$gene_id==gene & metab.new$hits > 0),]$name
bin.info.gene<-bin.info[which(bin.info$me %in% thebins),]
mean(bin.info.gene[phylum=="Actinobacteria"]$afe.median)
