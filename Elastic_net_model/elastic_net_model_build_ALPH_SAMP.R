#!/usr/bin/env Rscript


#############################################################
#Elastic net model building
#############################################################
#Requires
#R
#library glmnet
#Version 2.0.10 used for main manuscript results
#Newer versions will give very slightly different epigenetic age predictions (e.g. within a few hundredths). 

#This script builds leave-one-out models designed for parallelizing across samples (SAMP) and alpha values (ALPH).
#Each combination of sample and alpha value becomes it's own R script.
library(glmnet)

#############
#Read in data
#############
#Read in imputed mratios (output from impute.R)
epi<-read.table("./all_mratios_imputed_450k_n277_cleaned.txt",header=T,row.names=1)
head(epi)

#Read in metainfo with known chronological ages
info<-read.table("./info_for_age_predictions.txt",header=T)
head(info)

#length(which(colnames(epi2)==info2$Source.File.Name))

##########################
#Data Normalization
##########################
#Using alpha=ALPH
#Quantile normalize mratios by column (sample)
norm_counts<-as.matrix(apply(epi,2,function(x){return(qqnorm(x,plot=F)$x)}))

#Remove test subject(s)
#SAMP indexes from 1 to N samples
norm_train<-norm_counts[,-SAMP]
norm_test<-norm_counts[,SAMP]

#Quantile normalize training samples by row (CpG site)
trainreads_norm<-as.matrix(apply(norm_train,1,function(x){return(qqnorm(x,plot=F)$x)}))

#Create a vector of training and test ages for elastic net model construction
trainage<-info$age[-SAMP]
testage<-info$age[SAMP]

#QQ normalize each row in the test sample
#Note that this could be much more efficient using parallelR (e.g. parSapply)

#Store a new vector for normalizing the test sample
testreads_normalized<-norm_test

#For each CpG site
for (d in 1:length(testreads_normalized)){
  #Define the ECDF (depending on training sample size, this can be replaced with simply the training data -- norm_train)
  a<-ecdf(norm_counts[d,])
  #From this ECDF, return the probability of values being less than the training sample
  probs<-a(norm_test[d])
  #To avoid extreme values outside the range of the training samples, give 0's and 1's a manageable quantile (e.g. 0.99 and .01)
  #Note depending on the sample size, consider changing this number (e.g. to 1/N and 1-1/N respectively)
  probs[probs==1]<-.99
  probs[probs==0]<-.01
  #Given this probability, return the associated value from a standard normal that falls into the same quantile
  testreads_normalized[d]<-qnorm(probs)
}



###########################
#Elastic-net model building
###########################
#Using N-fold internal CV, train the elastic net model using the training data
#Note with larger sample sizes, N-fold internal CV becomes intractable
model<-cv.glmnet(trainreads_norm,trainage,nfolds=276,alpha=ALPH,standardize=F)

#Predict age using the test sample from parameters that minimized MSE during internal CV
predicted_SAMP<-predict(model,newx=t(testreads_normalized),s="lambda.min")

#Extract weights for this model
weights_SAMP<-unlist(coef(model,lambda="lambda.min"))[,1]

#Write out results for later concatenation
ifelse(!dir.exists(paths = "./predicted"),dir.create("./predicted"),FALSE)
ifelse(!dir.exists(paths = "./weights"),dir.create("./weights"),FALSE)

write.table(predicted_SAMP,"./predicted/predicted_QQ_n277_ALPH_SAMP.txt",quote=F,row.names=F,col.names=F)
write.table(weights_SAMP,"./weights/weights_QQ_n277_ALPH_SAMP.txt",quote=F,row.names=F,col.names=F)
