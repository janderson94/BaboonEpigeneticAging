#!/usr/bin/env Rscript

#Again this script is designed to run one chromosome (CHROM) at a time.
#This R script further filters the sites from get_data_CHROM.sh based on mean methratio and coverage.
#Requires data.table 
library(data.table)
data=fread('all_mratios_CHROMNAME_v2.txt',header=F)
data[,V9:=paste(data[,V1],data[,V2],sep="_")]

setkey(data,V9,V8)
temp=data
data=data.table(unique(data))
names=sort(unique(data[,V8]))
sites=unique(data[,V9])

#generate counts table
counts=as.data.table(unique(data[,V9]))
setnames(counts,"V1","site")
for(i in 1:length(names))
{
print(i)
counts[,names[i] := data[.(as.factor(counts[,site]),names[i]),V5]]
}

#generate mcounts table
mcounts=as.data.table(unique(data[,V9]))
setnames(mcounts,"V1","site")
for(i in 1:length(names))
{
print(i)
mcounts[,names[i] := data[.(as.factor(mcounts[,site]),names[i]),V6]]
}

counts[is.na(counts)] <- 0
mcounts[is.na(mcounts)] <- 0
mratios<-counts
mratios[,2:(length(names)+1)]<-mcounts[,2:(length(names)+1)]/counts[,2:(length(names)+1)]

#Calculate mean mratio, sd, and avg coverage
#Note that the SD filter here is not stringent, and could be much more strict.
#Depending on sample size and computational capacity, more stringent filtering could be done.
#Alternatively, if CpG sites are low quality (e.g. low depth, low variance), they will simply
mratios$avg=apply(mratios[,2:(length(names)+1)],1,mean,na.rm=TRUE)
mratios$sd=apply(mratios[,2:(length(names)+1)],1,sd,na.rm=TRUE)
mratios$depth<-apply(counts[,2:(length(names)+1)],1,mean,na.rm=TRUE)

#Filte
mratios2<-subset(mratios,avg<.9 & avg>.1 & sd > quantile(mratios$sd,probs=.05) & depth >5)
mcounts2<-mcounts[site %in%mratios2$site]
counts2<-counts[site %in%mratios2$site]

# save raw methylation ratios
write.table(mratios2,"methratio_CHROMNAME.txt",row.names=F,col.names=F,sep="\t")

# write methylated and total counts
write.table(mcounts2,file="mcounts_table_CHROMNAME.txt",row.names=F,sep="\t")
write.table(counts2,file="counts_table_CHROMNAME.txt",row.names=F,sep="\t")

# write mcounts and counts summaries, avg and sd
info=cbind(counts2$site,mratios2$avg,mratios2$sd,mratios2$depth)
# order = site ,mean, sd, depth
write.table(info,"info_CHROMNAME.txt",row.names=F,sep="\t")
