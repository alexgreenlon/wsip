library(data.table)

scaff.table <- fread("/Users/alexgreenlon/soil-C-SFA/wsip/results/dram/metawrap-drep-bins/flagella-scaffolds.csv")

assemblies <- unique(scaff.table$assembly)
assemblies <- assemblies[c(1:6,9:18)]

ass <- assemblies[1]
new.scaff.table <- scaff.table[assembly==ass]
ass.gc <- fread(paste("wsip/results/assembly-scaffold-stats/" ,ass ,"_scaffold_min1000.gc.txt",sep=""))
ass.gc[, c("scaffold","read","count") := tstrsplit(Name," ", fixed=TRUE)]
new.scaff.table <- merge(new.scaff.table,ass.gc,all.x=T)

for (i in 2:length(assemblies)){
  ass <- assemblies[i]
  ass.scaff.table <- scaff.table[assembly==ass]
  ass.gc <- fread(paste("wsip/results/assembly-scaffold-stats/" ,ass ,"_scaffold_min1000.gc.txt",sep=""))
  ass.gc[, c("scaffold","read","count") := tstrsplit(Name," ", fixed=TRUE)]
  ass.scaff.table <- merge(ass.scaff.table,ass.gc,all.x=T)
  new.scaff.table <- rbind(new.scaff.table,ass.scaff.table)
}

write.csv(new.scaff.table,"wsip/results/dram/metawrap-drep-bins/flagella-scaffold-lengths.csv",row.names=F,quote=F)

ass <- "H1-16-all-fractions"
new.new.scaff.table <- scaff.table[assembly==ass]
ass.gc <- fread(paste("wsip/results/assembly-scaffold-stats/" ,ass ,"_scaffold_min1000.gc.txt",sep=""))
ass.gc[, c("scaffold","flag","multi","len","read","count") := tstrsplit(Name," ", fixed=TRUE)]
new.new.scaff.table <- merge(new.new.scaff.table,ass.gc,all.x=T)

ass <- "H1-18-all-fractions"
ass.scaff.table <- scaff.table[assembly==ass]
ass.gc <- fread(paste("wsip/results/assembly-scaffold-stats/" ,ass ,"_scaffold_min1000.gc.txt",sep=""))
ass.gc[, c("scaffold","flag","multi","len","read","count") := tstrsplit(Name," ", fixed=TRUE)]
new.new.scaff.table <- merge(new.new.scaff.table,ass.gc,all.x=T)

write.csv(new.new.scaff.table,"wsip/results/dram/metawrap-drep-bins/flagella-scaffold-lengths.H1.csv",row.names=F,quote=F)

# manually combined the two sheets in google sheets, and then did some summary stats.
# https://docs.google.com/spreadsheets/d/1ESju1pf0W5a7pKGmtZKBf-GCEy7b_T5XGqddQsrFZQc/edit#gid=409818288
