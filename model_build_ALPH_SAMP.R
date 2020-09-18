#!/usr/bin/env Rscript
#SBATCH --get-user-env

#library(data.table)
library(glmnet)

epi<-read.table("all_methratios_imputed_450k_n277.txt",header=T,row.names=1)
colnames(epi)<-gsub(colnames(epi),pattern="^X",replacement="",fixed=F)

info<-read.table("info_n277_8jul19.txt",header=T)
epi<-epi[,colnames(epi) %in% info$file_name]

epi<-epi[,order(colnames(epi))]
info<-info[order(info$file_name),]

length(which(colnames(epi)==info$file_name))

#Using alpha=ALPHA
norm_counts<-as.matrix(apply(epi,2,function(x){return(qqnorm(x,plot=F)$x)}))

  #Remove test subject
  norm_train<-norm_counts[,-SAMP]
  norm_test<-norm_counts[,SAMP]

  trainreads_norm<-as.matrix(apply(norm_train,1,function(x){return(qqnorm(x,plot=F)$x)}))

  trainage<-info$age[-SAMP]
  testage<-info$age[SAMP]

  #QQ normalize each row
  testreads_normalized<-norm_test
  for (d in 1:length(testreads_normalized)){
    a<-ecdf(norm_counts[d,])
    probs<-a(norm_test[d])
    probs[probs==1]<-.99
    probs[probs==0]<-.01
    testreads_normalized[d]<-qnorm(probs)
  }

  model<-cv.glmnet(trainreads_norm,trainage,nfolds=276,alpha=ALPH,standardize=F)
  predicted_SAMP<-predict(model,newx=t(testreads_normalized),s="lambda.min")
  weights_SAMP<-unlist(coef(model,lambda="lambda.min"))[,1]

write.table(predicted_SAMP,"predicted_new/predicted_QQ_n277_ALPH_SAMP.txt",quote=F,row.names=F,col.names=F)
write.table(weights_SAMP,"weights_new/weights_QQ_n277_ALPH_SAMP.txt",quote=F,row.names=F,col.names=F)
