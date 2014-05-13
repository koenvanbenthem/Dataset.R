dat<-read.table(file="pop.csv",header=T)
dat2<-dat[,c("ID","t","z")]
library('reshape')
dat3<-reshape(dat2, idvar = "ID", timevar = "t", direction = "wide")
temp<-dat[,c(-1,-3,-5)]
temp<-temp[!duplicated(temp),]
dat4<-merge(dat3,temp,by="ID")
dat4