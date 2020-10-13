#!/usr/bin/env Rscript

##################################################################
#Filtering round 2
##################################################################
#Requires:
#R
#Data.table

#Again this script is designed to run one chromosome (CHROMNAME) at a time.
#This R script further filters the sites from get_data_CHROMNAME.sh based on mean methratio and coverage across samples.

###############################################
#Read in data output from get_data_CHROMNAME.sh
###############################################
library(data.table)
data=fread('all_mratios_CHROMNAME_v2.txt',header=F)
data[,V9:=paste(data[,V1],data[,V2],sep="_")]

setkey(data,V9,V8)
temp=data
data=data.table(unique(data))
names=sort(unique(data[,V8]))
sites=unique(data[,V9])

###################################
#Generate a counts table
###################################
#Rows are CpG sites
#Columns are samples
counts=as.data.table(unique(data[,V9]))
setnames(counts,"V1","site")
for(i in 1:length(names)){
counts[,names[i] := data[.(as.factor(counts[,site]),names[i]),V5]]
}

###################################
#Generate a methylated counts table
###################################
#Rows are CpG sites
#Columns are samples
mcounts=as.data.table(unique(data[,V9]))
setnames(mcounts,"V1","site")
for(i in 1:length(names)){
mcounts[,names[i] := data[.(as.factor(mcounts[,site]),names[i]),V6]]
}

#For sites with no data, give placeholder 0s for counts and methylated counts
counts[is.na(counts)] <- 0
mcounts[is.na(mcounts)] <- 0

###################################
#Generate a methratio table
###################################
mratios<-counts
#Methratios are methylated counts/total counts (i.e. mcounts/counts)
mratios[,2:(length(names)+1)]<-mcounts[,2:(length(names)+1)]/counts[,2:(length(names)+1)]

###################################
#Filter
###################################
#First, calculate mean mratio and avg coverage
mratios$avg=apply(mratios[,2:(length(names)+1)],1,mean,na.rm=TRUE)
mratios$depth<-apply(counts[,2:(length(names)+1)],1,mean,na.rm=TRUE)

#Note that the filters here are not stringent, and could be much more strict.
#Depending on sample size and computational capacity, more stringent filtering could be done.
#Alternatively, if CpG sites are low quality (e.g. low coverage), they will simply not be weighted in the model if using for elastic net modeling
#Get per-CpG site average methratio
#Get per-CpG site average depth

#Filter for mean methratio between 0.1 and 0.9 and average depth >5
mratios2<-subset(mratios,avg<.9 & avg>.1 & depth >5)
mcounts2<-mcounts[site %in%mratios2$site]
counts2<-counts[site %in%mratios2$site]

###################################
# Save raw methylation ratios
###################################
write.table(mratios2,"methratio_CHROMNAME.txt",row.names=F,col.names=F,sep="\t")

###################################
# write methylated and total counts
###################################
write.table(mcounts2,file="mcounts_table_CHROMNAME.txt",row.names=F,sep="\t")
write.table(counts2,file="counts_table_CHROMNAME.txt",row.names=F,sep="\t")

# write out summaries if desired
info=cbind(counts2$site,mratios2$avg,mratios2$depth)
# order = site ,mean mratio, and depth
write.table(info,"info_CHROMNAME.txt",row.names=F,sep="\t")
