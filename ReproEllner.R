setwd(dir="C:/Users/Timothée/Documents/GitHub/Dataset.R/Data/simple")
popfile<-read.table(file="TheDataSetV0.csv",header=T)
popfile$E<-popfile$z-popfile$bvs

meanZ<-tapply(popfile$z[which(popfile$phi==1)],INDEX = popfile$t[which(popfile$phi==1)],FUN = mean)
meanBVS<-tapply(popfile$bvs[which(popfile$phi==1)],INDEX = popfile$t[which(popfile$phi==1)],FUN = mean)
lmZ<-lm(meanZ~1+I(1:length(meanZ)))
lmBVS<-lm(meanBVS~1+I(1:length(meanBVS)))
meanE<-meanZ-meanBVS

lmE<-lm(meanE~1+I(1:length(meanE)))

meanARS<-tapply(popfile$ARS[which(popfile$phi==1)],INDEX = popfile$t[which(popfile$phi==1)],FUN = mean)
lmARS<-summary(glm(meanARS~1+meanE+meanBVS))

E1<-coef(lmE)[1]+coef(lmE)[2]*1
E50<-coef(lmE)[1]+coef(lmE)[2]*50

BVS1<-coef(lmBVS)[1]+coef(lmBVS)[2]*1
BVS50<-coef(lmBVS)[1]+coef(lmBVS)[2]*50

S_E1_BVS1<-E1*coef(lmARS)["meanE",1]+BVS1*coef(lmARS)["meanBVS",1]+coef(lmARS)["(Intercept)",1]
S_E1_BVS50<-E1*coef(lmARS)["meanE",1]+BVS50*coef(lmARS)["meanBVS",1]+coef(lmARS)["(Intercept)",1]
S_E50_BVS1<-E50*coef(lmARS)["meanE",1]+BVS1*coef(lmARS)["meanBVS",1]+coef(lmARS)["(Intercept)",1]
S_E50_BVS50<-E50*coef(lmARS)["meanE",1]+BVS50*coef(lmARS)["meanBVS",1]+coef(lmARS)["(Intercept)",1]

evol<-0.5*(S_E1_BVS50-S_E1_BVS1)+0.5*(S_E50_BVS50-S_E50_BVS1)
ecol<-0.5*(S_E50_BVS1-S_E1_BVS1)+0.5*(S_E50_BVS50-S_E1_BVS50)


#year by year
Farse<-function(popfile,year){
  meanZ<-tapply(popfile$z[which(popfile$phi==1)],INDEX = popfile$t[which(popfile$phi==1)],FUN = mean)
  meanBVS<-tapply(popfile$bvs[which(popfile$phi==1)],INDEX = popfile$t[which(popfile$phi==1)],FUN = mean)
  lmZ<-lm(meanZ~1+I(1:length(meanZ)))
  lmBVS<-lm(meanBVS~1+I(1:length(meanBVS)))
  meanARS<-tapply(popfile$ARS[which(popfile$phi==1)],INDEX = popfile$t[which(popfile$phi==1)],FUN = mean)
  lmARS<-summary(glm(meanARS~1+meanZ+meanBVS))
  Z1<-coef(lmZ)[1]+coef(lmZ)[2]*1
  Z30<-coef(lmZ)[1]+coef(lmZ)[2]*30
  BVS1<-coef(lmBVS)[1]+coef(lmBVS)[2]*1
  BVS30<-coef(lmBVS)[1]+coef(lmBVS)[2]*30
  S_Z1_BVS1<-Z1*coef(lmARS)["meanZ",1]+BVS1*coef(lmARS)["meanBVS",1]+coef(lmARS)["(Intercept)",1]
  S_Z1_BVS30<-Z1*coef(lmARS)["meanZ",1]+BVS30*coef(lmARS)["meanBVS",1]+coef(lmARS)["(Intercept)",1]
  S_Z30_BVS1<-Z30*coef(lmARS)["meanZ",1]+BVS1*coef(lmARS)["meanBVS",1]+coef(lmARS)["(Intercept)",1]
  S_Z30_BVS30<-Z30*coef(lmARS)["meanZ",1]+BVS30*coef(lmARS)["meanBVS",1]+coef(lmARS)["(Intercept)",1]
  evol<-0.5*(S_Z1_BVS30-S_Z1_BVS1)+0.5*(S_Z30_BVS30-S_Z30_BVS1)
  ecol<-0.5*(S_Z30_BVS1-S_Z1_BVS1)+0.5*(S_Z30_BVS30-S_Z1_BVS30)
}
