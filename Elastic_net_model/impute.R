#!/usr/bin/env Rscript

#############################################################
#Impute missing methratio calls prior to elastic net modeling
#############################################################
#Requires 
#R
#libraries impute and data.table
library(data.table)
library(impute)

#Load in the combined set of sites output from get_data_CHROM.R across, concatenated into one file "all_methratios.txt"
all<-fread("all_methratios.txt")



##################################################
#Final filter based on missingness for elastic net
##################################################

#Calculate the number of missing methratio calls per site
missing<-apply(all,1,function(x){length(which(is.na(x)))})

#Remove sites missing in >5% of the data
all2<-all[missing<14,]
rownames(all2)<-all2$site; all2$site<-NULL



##################################################
#Impute using default k-nearest neighbors approach
##################################################
all_imputed<-impute.knn(as.matrix(all2))$data

#Write out impute mratios for elastic net modeling
write.table(all_imputed,"all_methratios_imputed.txt",quote=F,row.names=T,col.names=T)

