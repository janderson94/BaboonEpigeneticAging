#!/usr/bin/env Rscript

###################################################################################
################################## Main results ###################################
###################################################################################

#################################################
#Read in meta info and epigenetic age predictions
#################################################
#This info file is SI table 1 with renamed columns for readability
#N277 model
info<-read.csv("./SI_table1.csv")[1:277,1:14]


############################
#How well do we predict age?
############################
#In terms of median absolute difference
median(abs(info$Predicted_epi_age-info$Age))
sd(abs(info$Predicted_epi_age-info$Age))

#In terms of pearson's correlation
#Line 97 is rounded to 9.71 x 10 ^(-54)
out<-cor.test(info$Predicted_epi_age,info$Age);out$estimate;out$p.value


###############
#Repeat samples
###############
#Id indexes individuals
info$id<-gsub(gsub(info$SRA.sample.ID,pattern="_B",replacement="",fixed=T),pattern="_C",replacement = "",fixed=T)

#30individuals have repeat samples
length(table(info$id)[table(info$id)>1])
info2<-info[info$id %in% names(table(info$id)[table(info$id)>1]),]
info2<-info2[order(info2$SRA.sample.ID,info2$Age),]

#For two individuals with 3 samples in the dataset, drop the middle sample for this analysis.
info2<-info2[-which(rownames(info2)%in% c('86','149')),]

#In 26/30 sample pairs, the older sample pair is predicted to be older than the younger sample pair
length(which(info2$Predicted_epi_age[seq(2,60,2)]>info2$Predicted_epi_age[seq(1,59,2)]))

binom.test(x = 26,n = 30)

#Does the difference in age between resamplings predict the difference in predicted age, controlling for sex?
dif_pred<-info2$Predicted_epi_age[seq(2,60,2)]-info2$Predicted_epi_age[seq(1,59,2)]
dif_age<-info2$Age[seq(2,60,2)]-info2$Age[seq(1,59,2)]
sex<-info2$Sex[seq(2,60,2)]

summary(lm(dif_pred~dif_age+sex))

mean(dif_pred-dif_age)
sd(dif_pred-dif_age)



#####################################
#Sex differences in clock performance
#####################################

####Males
#Median absolute difference
m<-which(info$Sex=="M")
median(abs(info$Predicted_epi_age[m]-info$Age[m]))
sd(abs(info$Predicted_epi_age[m]-info$Age[m]))

#Pearson's correlation
out<-cor.test(info$Predicted_epi_age[m],info$Age[m]);out$estimate;out$p.value


#Females
#Median absolute difference
f<-which(info$Sex=="F")
median(abs(info$Predicted_epi_age[f]-info$Age[f]))
sd(abs(info$Predicted_epi_age[f]-info$Age[f]))

#Pearson's correlation
out<-cor.test(info$Predicted_epi_age[f],info$Age[f]);out$estimate;out$p.value


#Sex differences in absolute differences
wilcox.test(x = abs(info$Predicted_epi_age[m]-info$Age[m]),
            abs(info$Predicted_epi_age[f]-info$Age[f]))


#Is there a significant difference in terms of slope of epigenetic age trajectories?
summary(lm(Predicted_epi_age~Age*Sex,data=info))->b;b$coefficients


#Is this sex difference apparent in younger individuals?
summary(lm(Predicted_epi_age~Age*Sex,data=info,subset=which(info$Age<8)))->out;out$coefficients

#No, only after they have reached physiological and social adulthood
summary(lm(Predicted_epi_age~Age*Sex,data=info,subset=which(info$Age>8)))->out;out$coefficients

rm(out)

######################################################################
#What predicts deviations in age predictions?
######################################################################
#Here we explore the ability of 4 fitness-predictive factors to predict epigenetic age deviations
#These include:
#BMI - Body Mass Index
#Cumulative Early Adversity Score
#Dyadic Social Connectedness (Female to Female bond strengths)
#Ordinal dominance rank

#Our outcome variable is simply the difference between predicted epigenetic age and known chronological age
info$delta_age<-info$Predicted_epi_age-info$Age


#First we note that dominance rank in males and BMI are correlated
cor.test(info$BMI[info$Sex=="M"],
         info$Rank[info$Sex=="M"],use = "pairwise.complete.obs")


#Build piecewise regression to generate age-adjusted BMI
#Requires segmented and scales packages
library(segmented)
library(scales)

info$Age.adjusted.BMI_new<-NA
info_males<-info[m,]

fit.lm.male<-glm(BMI~Age,family="gaussian",data = info_males)

fit.seg.male<-segmented(fit.lm.male,seg.Z=~Age,psi=8)
predicted<-predict(object=fit.seg.male,newdata = info_males)

par(mfrow=c(1,1),pty="s")
#Plot our segmented fit line
plot.segmented(fit.seg.male,ylim=c(25,60),xlim=c(0,20),xlab="",ylab="")
par(new=T)
#Overlay our actual data
plot(info_males$BMI~info_males$Age,subset = m,col=alpha("Steel Blue",alpha=0.8),pch=20,
     xlab="Age (years)",ylab="Body Mass Index (Kg/m^2)",main="Males",cex=1.5,ylim=c(25,60),xlim=c(0,20))
par(new=T)
#Finally plot our predicted BMI (based on ages) as points
plot(predicted~info_males$Age,ylim=c(25,60),xlim=c(0,20),xlab="",ylab="")

#Age adjusted BMI is simply the residuals between BMI from the best-fit piecewise model
info_males$Age.adjusted.BMI_new<-info_males$BMI-predicted

#The correlations between the newly generated age-adjusted BMI is effectively 1. 
#But there is slight stochasticity in these results
plot(info_males$Age.adjusted.BMI_new~info_males$Age.adjusted.BMI)
summary(lm(info_males$Age.adjusted.BMI_new~info_males$Age.adjusted.BMI))
mean(abs(info_males$Age.adjusted.BMI_new-info_males$Age.adjusted.BMI),na.rm=T)

#We suggest using our stored age-adjusted BMI values in info to recreate our exact results



############################################
#Cross-sectional analyses
############################################

#Age-adjusted BMI breaks the correlation between Male BMI and rank
cor.test(info$Age.adjusted.BMI[info$Sex=="M"],
    info$Rank[info$Sex=="M"],use = "pairwise.complete.obs")



########################################
################ FEMALES ###############
########################################

summary(lm(info$delta_age~info$Cumulative.Early.Adversity.Score+
             info$Dyadic.Social.Connectedness..F.to.F.+
             info$Rank+
             info$Age.adjusted.BMI+
             info$Age,subset = f))

#Note for p-values < 2e-16, the exact value can be recalled via:
out<-summary(lm(info$delta_age~info$Cumulative.Early.Adversity.Score+
             info$Dyadic.Social.Connectedness..F.to.F.+
             info$Rank+
             info$Age.adjusted.BMI+
             info$Age,subset = f))
out$coefficients
rm(out)


########################################
################# MALES ################
########################################

summary(lm(info$delta_age~info$Cumulative.Early.Adversity.Score+
             info$Rank+
             info$Age.adjusted.BMI+
             info$Age,subset = m))


###################################
#Alternate male models, SI Table 5
###################################
info_males<-info[m,]

#Modeling rank after controling for BMI
info_males$rank_no_bmi<-NA

info_males$rank_no_bmi[as.numeric(names(residuals(lm(info_males$Rank~info_males$BMI))))]<-residuals(lm(info_males$Rank~info_males$BMI))

summary(lm(info_males$delta_age~info_males$Cumulative.Early.Adversity.Score+
             info_males$rank_no_bmi+
             info_males$Age.adjusted.BMI+
             info_males$Age))

#What if we modeled raw BMI (i.e. not age-adjusted BMI)?
summary(lm(info_males$delta_age~info_males$Cumulative.Early.Adversity.Score+
             info_males$Rank+
             info_males$BMI+
             info_males$Age))

#Or if we modeled BMI controlling for rank?
info_males$bmi_no_rank<-NA
info_males$bmi_no_rank[as.numeric(names(residuals(lm(info_males$BMI~info_males$Rank))))]<-residuals(lm(info_males$BMI~info_males$Rank))

summary(lm(info_males$delta_age~info_males$Cumulative.Early.Adversity.Score+
             info_males$Rank+
             info_males$bmi_no_rank+
             info_males$Age))

#What if we exclude small (phyiscally immature) males?
large<-which(info_males$BMI>41)

summary(lm(info_males$delta_age~info_males$Cumulative.Early.Adversity.Score+
             info_males$Rank+
             info_males$Age.adjusted.BMI+
             info_males$Age,subset = large))

#Doing so breaks the raw BMI/rank correlation but still finds a significant rank effect
cor.test(info_males$Rank[large],info_males$BMI[large],
    use="pairwise.complete.obs")


#What if we control for the non-linear age structure of rank?
#Is rank-for-age a better predictor than absolute rank?
info_males$age_squared<-(info_males$Age^2)
info_males$rank_no_age<-NA
info_males$rank_no_age[as.numeric(names(residuals(lm(info_males$Rank~info_males$Age))))]<-
  residuals(lm(info_males$Rank~info_males$Age+ info_males$age_squared))

#No-- rank remains marginally significant, and rank-for-age is not found to be significant
summary(lm(info_males$delta_age~info_males$Cumulative.Early.Adversity.Score+
             info_males$Rank+
             info_males$Age.adjusted.BMI+
             info_males$Age+
             info_males$rank_no_age))



############################################
#Longitudinal analyses
############################################
#Here we generated additional data from males already in the dataset to better ask about longitudinal changes in dominance rank and epigenetic age
info<-read.csv("./SI_table1.csv")

#Becasue of the cross-sectional results, we are focusing on males here.
info<-info[info$Sex=="M",]

#Calculating rank_for_age for later
info$age_squared<-(info$Age^2)
info$rank_no_age<-NA
info$rank_no_age[as.numeric(names(residuals(lm(info$Rank~info$Age))))]<-
  residuals(lm(info$Rank~info$Age+info$age_squared))

#Calculating rank controlling for BMI for later
info$rank_no_bmi<-NA
info$rank_no_bmi[as.numeric(names(residuals(lm(info$Rank~info$BMI))))]<-residuals(lm(info$Rank~info$BMI))

#Calculating BMI controlling for rank for later
info$bmi_no_rank<-NA
info$bmi_no_rank[as.numeric(names(residuals(lm(info$BMI~info$Rank))))]<-residuals(lm(info$BMI~info$Rank))


#For these samples, we will use the age predictions from the N286 model. 
info$delta_age<-info$Predicted_epi_age_n286-info$Age

#Across all males, regress out the age bias in delta age
info$residual_age<-residuals(lm(info$delta_age~info$Age))

#Again, ID indexes the individuals
info$id<-gsub(gsub(info$SRA.sample.ID,pattern="_B",replacement="",fixed=T),pattern="_C",replacement = "",fixed=T)

#Exclude males without rank values
info<-info[!is.na(info$Rank),]

#Retain individuals with multiple samples
info2<-info[info$id %in% names(table(info$id)[table(info$id)>1]),]

#Order repeat samples by age at sampling
info2<-info2[order(info2$SRA.sample.ID,info2$Age),]

#For the four males with three samples in the dataset, drop the middle sample, data permitting
table(info2$id)[table(info2$id)>2]

#Which samples to drop
#AMB_282 does not have rank associated with one of the paired samples
#Exclude AMB 69's middle sample 
#Exclude AMB_133's middle sample
#Exclude AMB_152's middle sample
#Exclude the middle AMB 230 sample (all rank 1)

info2<-info2[-which(info2$SRA.sample.ID %in% c('AMB_282','AMB_282_B',
                                            'AMB_230_B',
                                            "AMB_133_B",
                                            "AMB_152_B",
                                            "AMB_69_C")),]

#This leaves us with 14 males, each with 2 samples.
#13/14 males occupied different ranks between the two samplings
length(unique(info2$id))
info3<-info2[order(info2$id,info2$Rank),]

#Exclude the male who was rank 1 for both sample collections
info4<-info3[info3$id != "AMB_230",]

#Are residual epigenetic ages higher when the samples were taken from males when they occupied higher ranks?
t.test(info4$residual_age[seq(1,25,2)],info4$residual_age[seq(2,26,2)],paired=T)
rm(info4)


#Does change in rank significantly predict change in residual epigenetic age?
#Here we can use the male who was rank 1 from both sample collections
info3<-info2[order(info2$id,info2$Age),]

dif_rank<-info3$Rank[seq(2,28,2)]-info3$Rank[seq(1,27,2)]
dif_residual_age<-info3$residual_age[seq(2,28,2)]-info3$residual_age[seq(1,27,2)]

#Yes, it does.
summary(lm(dif_residual_age~dif_rank))


#Differences in BMI
dif_bmi<-info3$BMI[seq(2,28,2)]-info3$BMI[seq(1,27,2)]
dif_residual_age<-info3$residual_age[seq(2,28,2)]-info3$residual_age[seq(1,27,2)]

summary(lm(dif_residual_age~dif_bmi))


#Rank controlling for BMI
dif_rank_no_bmi<-info3$rank_no_bmi[seq(2,28,2)]-info3$rank_no_bmi[seq(1,27,2)]

summary(lm(dif_residual_age~dif_rank_no_bmi))


#BMI no rank
dif_bmi_no_rank<-info3$bmi_no_rank[seq(2,28,2)]-info3$bmi_no_rank[seq(1,27,2)]

summary(lm(dif_residual_age~dif_bmi_no_rank))


##################################################################
###################### End main results ##########################
##################################################################
