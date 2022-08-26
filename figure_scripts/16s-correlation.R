dat<-read.csv("soil-C-SFA/wsip/results/16S-comparison/all.16s-all-contigs.bs.test.csv")
head(dat)
lm.contig<-lm(X16S.afe.median ~ AFE, data=dat)
summary(lm.contig)
head(dat)
lm.bin<-lm(X16S.afe.median ~ bin.afe.median, data=dat)
summary(lm.bin)
history(p<0.05
history()
