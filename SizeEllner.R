setwd(dir="C:/Users/Timothée/Documents/GitHub/Dataset.R/Data/simple")
popfile<-read.table(file="TheDataSetV0.csv",header=T)

popfile$birthSize<-NA
IDS<-unique(popfile$ID)
for (i in 1:length(IDS))
{
  popfile$birthSize[which(popfile$ID==IDS[i])]<-popfile$z[which(popfile$ID==IDS[i] & (popfile$age==0 | popfile$t==1))]
}

popfile$Growth<-NA
for (i in 1:nrow(popfile))
  {
    Zt<-popfile$z[i]
    time<-popfile$t[i]-1

    if(nrow(popfile[which(popfile$ID==popfile$ID[i] & popfile$t==time),])>0)
      {
        Ztt<-popfile$z[which(popfile$ID==popfile$ID[i] & (popfile$t)==time)]
        popfile$Growth[i]<-Zt-Ztt
      }
  }

meanG<-tapply(popfile$Growth[which(popfile$phi==1)],INDEX = popfile$t[which(popfile$phi==1)],FUN = function(x){mean(x,na.rm=T)})
meanP<-tapply(popfile$birthSize[which(popfile$phi==1)],INDEX = popfile$t[which(popfile$phi==1)],FUN = function(x){mean(x,na.rm=T)})
meanBVS<-tapply(popfile$bvs[which(popfile$phi==1)],INDEX = popfile$t[which(popfile$phi==1)],FUN = mean)
meanE<-meanP-meanBVS

lmP<-lm(meanP~1+I(1:length(meanP)))
lmG<-lm(meanG~1+I(1:length(meanG)))
lmBVS<-lm(meanBVS~1+I(1:length(meanBVS)))
lmE<-lm(meanE~1+I(1:length(meanE)))

meanZ<-tapply(popfile$z[which(popfile$phi==1)],INDEX = popfile$t[which(popfile$phi==1)],FUN = mean)
lmZ<-summary(glm(meanZ~1+meanE+meanBVS+meanG))
meanE[-1]+meanG[-1]+meanBVS[-1]-meanZ[-1]
summary(glm(meanZ~1+meanBVS))

E1<-coef(lmE)[1]+coef(lmE)[2]*1
E50<-coef(lmE)[1]+coef(lmE)[2]*50

BVS1<-coef(lmBVS)[1]+coef(lmBVS)[2]*1
BVS50<-coef(lmBVS)[1]+coef(lmBVS)[2]*50

S_E1_BVS1<-E1*coef(lmZ)["meanE",1]+BVS1*coef(lmZ)["meanBVS",1]+coef(lmZ)["(Intercept)",1]
S_E1_BVS50<-E1*coef(lmZ)["meanE",1]+BVS50*coef(lmZ)["meanBVS",1]+coef(lmZ)["(Intercept)",1]
S_E50_BVS1<-E50*coef(lmZ)["meanE",1]+BVS1*coef(lmZ)["meanBVS",1]+coef(lmZ)["(Intercept)",1]
S_E50_BVS50<-E50*coef(lmZ)["meanE",1]+BVS50*coef(lmZ)["meanBVS",1]+coef(lmZ)["(Intercept)",1]

evol<-0.5*(S_E1_BVS50-S_E1_BVS1)+0.5*(S_E50_BVS50-S_E50_BVS1)
ecol<-0.5*(S_E50_BVS1-S_E1_BVS1)+0.5*(S_E50_BVS50-S_E1_BVS50)


rm(list=ls())
setwd("/home/koen/Dropbox/Decomposing pop dynamics/Koen/Coulson_Tulja_Sim")
dat<-read.table(file="TheDataSetV0.csv",header=T)
dat<-subset(dat,phi==1)
frame<-aggregate(dat$z,by=list(dat$t),mean)
colnames(frame)<-c('t','z')
frame$a<-aggregate(dat$bvs,by=list(dat$t),mean)[,2]

whereis<-function(x){
  if(x %in% dat$ID){
    return(min(which(dat$ID==x)))
  }else{
    return(NA)
  }
}
indices<-sapply(1:max(dat$ID),whereis)
dat$bw<-dat$z[indices[dat$ID]]
frame$bw<-aggregate(dat$bw,by=list(dat$t),mean)[,2]
frame$e<-frame$bw-frame$a
model_bw<-lm(bw~t,data=frame)
model_a<-lm(a~t,data=frame)
a<-function(t){
  predict(model_a,newdata=list(t=t),type="response")
}

e<-function(t){
  predict(model_bw,newdata=list(t=t),type="response")-a(t)
}

# Growth since last timestep
dat$GR<-NA
for(i in 1:length(dat$t)){
  indi<-(dat$ID==dat$ID[i] & dat$t==(dat$t[i]-1))
  if(sum(indi)==1){
    dat$GR[i]<-dat$z[i]-dat$z[indi]
  }else if(sum(indi)>1){
    warning("oops")
  }
}
frame$GR<-aggregate(dat$GR,by=list(dat$t),mean,na.rm=TRUE)[,2]
model_c<-lm(GR~t,data=frame)
model_z<-lm(z~a+e+GR,data=frame)
k<-function(t){
  predict(model_c,newdata=list(t=t),type="response")
}
z<-function(a,e,c){
  predict(model_z,newdata=list(a=a,e=e,GR=c),type="response")
}

framepje<-expand.grid(ai=c(0,1),ei=c(0,1),ci=c(0,1))
framepje$a<-a(49*framepje$ai+1)
framepje$e<-e(49*framepje$ei+1)
framepje$k<-k(49*framepje$ci+1)
framepje$z<-z(framepje$a,framepje$e,framepje$k)
lm(z~ai+ei+ci,data=framepje)