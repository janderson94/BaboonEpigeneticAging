#!/usr/bin/env Rscript

#############################################################
#Impute missing methratio calls prior to elastic net modeling
#############################################################
#Requires 
#R
#libraries impute and data.table

library(data.table)

all<-fread("all_methratios.txt")

missing<-apply(all,1,function(x){length(which(is.na(x)))})

#Want to remove sites missing in >5% of the data for the purposes of elastic net efficiency
all2<-all[missing<14,]
rownames(all2)<-all2$site
all2$site<-NULL

library(impute)

#Save the imputed, missing methylation ratios
all_imputed<-impute.knn(as.matrix(all2))$data

#Write out impute mratios for elastic net modeling
write.table(all_imputed,"all_methratios_imputed_450k_n277.txt",quote=F,row.names=T,col.names=T)

