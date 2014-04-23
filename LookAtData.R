setwd(dir="/GitHub/Dataset.R")
popfile<-read.table(file="pop.csv",header=T)
plot(popfile$z,x=popfile$age)
plot(popfile$z,x=popfile$C)
plot(popfile$,x=popfile$C)

library(lme4)
mm0<-lmer(z~1+age+t+s+(1|ID),data=popfile)
summary(mm0)
table(popfile$C)

Camemberts<-tapply(X=popfile$C,INDEX=popfile$t,FUN=sum)#retrieve the yearly food abundance
popSize<-tapply(X=popfile$ID[which(popfile$phi==1)],INDEX=popfile$t[which(popfile$phi==1)],FUN=length)#retrieve the yearly pop size

plot(Camemberts,popSize)

mPhiC0<-glm(phi~1+C+age+z,data=popfile)
summary(mPhiC0)
