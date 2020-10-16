####################################
#Simulations for reviewers
####################################

#This script documents and utilizes code to simulate a random variable with a prticular correlation to a vector of interest
#We use these simulations to show that the bias in age predictions must be 
#accounted for when testing for predictors of epigenetic age acceleration.

#Simcor is main function that takes a variable, x, and simulates a new, random variable with a specified:
#Mean, sd, and correlation to x

simcor <- function (x, ymean=MEAN, ysd=SD, correlation=0) {
  #What a vector of the same length
  n <- length(x)
  
  #Draw from a normal distribution based on specified mean and SD
  y <- rnorm(n,mean = ymean,sd=ysd)
  
  #Take the vector of interest, scaled to mean 0, var 1, rescaled by the correlation
  #Add Randomly distributed noise scaled by 1-correlation
  #Result will be mean 0, var 1 vector
  z <- correlation * scale(x)[,1] + sqrt(1 - correlation^2) * scale(resid(lm(y ~ x)))[,1]
  
  #Rescale by desired sd and add to desired mean
  yresult <- ymean + ysd * z
  return(yresult)
}


###############################################################
##The data frame "data" loaded in here is a slice of SI Table 1
##It focuses on males from the main N=277 model 
#Load the RData from the directory where this Simulation folder is cloned.
load("./.RData")
head(data)

#Here we show that if a variable is correlated with chronological age, like raw BMI is,
#then age bias in delta age must be account for.
#For example,
par(mfrow=c(1,1),pty="s",cex=1.5)
plot(data$delta_age~data$age,pch=20,xlab="Chronological age (years)",ylab="Delta age")
abline(lm(data$delta_age~data$age),lty=2)

#Here we see the strong age bias in our delta age predictions, such that older individuals tend to be underestimated
#i.e. they tend to have a negative delta age (Predicted-chronological age)

#Since BMI has a strong relationship with age throughout development
plot(data$BMI~data$age,pch=20,xlab="Chronological age (years)",ylab="BMI")

#The age structure in both delta age and BMI produce a correlation
plot(data$delta_age~data$BMI,pch=20,xlab="BMI",ylab="Delta age") 

############################################
#To show that this is a statistical artifact,
#and not a product of a true relationship between BMI and epigenetic age acceleration:

#We will simulate a random variable 1000 times that has the same mean and variance as BMI, 
#and this variable has the same correlation with age that BMI has with age 

#If these simulated variables predict delta age, despite having no true relationship, this result would suggest the association
#between BMI and delta age is also an artifact of it's correlation with chronological age

#If this simulated variable does not predict delta age, this result would suggest that 
#the association between BMI and delta age is due to a true biological relationship between variance in BMI and epigenetic age

#The correlation between BMI and age in these samples is ~0.6
cor(data$BMI,data$age,use="pairwise.complete.obs")

#Lets simulate one, random variable as outlined above.
set.seed(10)
data$simulated<-simcor(x = data$age,ymean=mean(data$BMI,na.rm = T),ysd=sd(data$BMI,na.rm=T),correlation = 0.60)

par(mfrow=c(1,3),pty="s")
plot(data$age~data$simulated,xlab="Simulated variable",ylab="Chronological age",pch=20)
cor.test(data$simulated,data$age)
#As expected, the simulated variable and age have a correlation of 0.6

#We also see that the simulated variable and delta age are strongly, negatively correlated (as is BMI and delta age)
plot(data$delta_age~data$simulated,xlab="Simulated variable",ylab="Delta age",pch=20)
cor.test(data$simulated,data$delta_age)

#Finally, using residual epigenetic age (delta age controlling for the age bias in predictions by regressing age out from delta age and taking the residuals),
#We find no relationship between the simulated variable and residual epigenetic age
plot(data$Residual_epi_age~data$simulated,xlab="Simulated variable",ylab="Residual epigenetic age",pch=20)
cor.test(data$Residual_epi_age,data$simulated)

#To make sure this is not a fluke, lets rerun this 1000 times.
#Each time we will simulate a new, random variable.
#We will extract the pearson's r correlation between the simulated variable and delta age
#We will extract the pearson's r correlation between the simulated variable and residual epigenetic age
#And we will extract the associated significance in those correlations as well.

#Simulate the variable 1000 times and take correlation with predicted-epi and residual epi age
b1<-1:1000
b2<-1:1000
p1<-1:1000
p2<-1:1000

for(f in 1:1000){
  s<-simcor(x = data$age,ymean=mean(data$BMI,na.rm = T),ysd=sd(data$BMI,na.rm=T),correlation = 0.60)
  b1[f]<-cor(s,data$delta_age,use="pairwise.complete.obs")
  b2[f]<-cor(s,data$Residual_epi_age,,use="pairwise.complete.obs")
  p1[f]<-unlist(cor.test(s,data$delta_age,use="pairwise.complete.obs")[3])
  p2[f]<-unlist(cor.test(s,data$Residual_epi_age,use="pairwise.complete.obs")[3])
  
}

#Visualize results
#Colored by significance
#First plot a histogram of the correlation between the simulated variable and delta age (i.e. predicted- chronological age)
library(scales)
par(mfrow=c(1,2),pty="s")
hist(b1[p1>.1],xlab="Pearson's r",breaks=25,main="Simulated variable vs Delta age", 
     xlim=c(min(b1),max(b1)),col=alpha("Black",alpha=.6),ylim=c(0,50)) 
par(new=TRUE)
hist(b1[p1<.05],xlab="Pearson's r",breaks=50,main="",
     xlim=c(min(b1),max(b1)),col=alpha("Steel Blue",alpha=.6),ylim=c(0,50)) 
par(new=TRUE)
hist(b1[p1<.1 & p1>.05],xlab="Pearson's r",breaks=5,main="",
     xlim=c(min(b1),max(b1)),col=alpha("#70869E",alpha=.6),ylim=c(0,50)) 
#Dashed line indicates the strength of the observed correlation between BMI and delta age
abline(v=cor(data$BMI,data$delta_age,use="pairwise.complete.obs"),lty=2,lwd=3)

#Many of the 1,000 simulated variables are significantly associated with delta age at a nominal p-value of 0.05
length(which(p1<.05))


#Plot the same results but from when we regress the age structure out of delta age first.
hist(b2[p2>.1],xlab="Pearson's r",breaks=25,main="Simulated variable vs Residual epi age",
     xlim=c(min(b2),max(b2)),col=alpha("Black",alpha=.6),ylim=c(0,80)) 
par(new=TRUE)
hist(b2[p2<.05],xlab="Pearson's r",breaks=50,main="",
     xlim=c(min(b2),max(b2)),col=alpha("Steel Blue",alpha=.6),ylim=c(0,80))
par(new=TRUE)
hist(b2[p2>.05 & p2<.1],xlab="Pearson's r",breaks=25,main="",
     xlim=c(min(b2),max(b2)),col=alpha("#70869E",alpha=.6),ylim=c(0,80))
#Dashed line indicates the strength of the observed correlation between BMI and residual epi age
abline(v=cor(data$BMI,data$Residual_epi_age,use="pairwise.complete.obs"),lwd=3,lty=2)

#Far fewer of the 1,000 simulated variables are significantly associated with delta age at a nominal p-value of 0.05
#(Closer to the 5% expectation)
length(which(p2<.05))

#rm(b1,b2,f,p1,p2,s)

my_session
#my_session<-sessionInfo()







