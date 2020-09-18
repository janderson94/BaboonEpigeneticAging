#!/usr/bin/env Rscript
#SBATCH --get-user-env

library(data.table)

all<-fread("all_methratios.txt")

temp<-read.table("all_methratios_no_missing_n277.txt",header=T)

all[,283:285]<-NULL

colnames(all)<-colnames(temp)

missing<-apply(all,1,function(x){length(which(is.na(x)))})

all2<-all[missing<14,]
rownames(all2)<-all2$site
all2$site<-NULL

library(impute,lib.loc="/data/tunglab/jaa57/programs/R_packages")

all_imputed<-impute.knn(as.matrix(all2))$data

dim(all_imputed)
write.table(all_imputed,"all_methratios_imputed_450k_n277.txt",quote=F,row.names=T,col.names=T)

