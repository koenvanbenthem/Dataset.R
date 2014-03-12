setwd(dir="/GitHub/Dataset.R")
popfile<-read.table(file="pop.csv",header=T)
plot(popfile$z,x=popfile$age)

library(lme4)
mm0<-lmer(z~1+age+t+s+(1|ID),data=popfile)
summary(mm0)
