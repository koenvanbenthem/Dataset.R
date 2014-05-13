setwd(dir="/GitHub/Dataset.R")
popfile<-read.table(file="pop.csv",header=T)
plot(popfile$z,x=popfile$age)

plot(popfile$z,x=popfile$C)#size depending on how many Camemberts you (leprechaun) get, or the other way around rather

library(ggplot2)
qplot(x=bvs,y=z,colour=t,data=popfile)#visualize heritability and response to selection (temporal change in breeding values)
plot(popfile$z,x=popfile$t)
plot(popfile$bvs,x=popfile$t)#there is not much genetic variance left after 15 years!

popsize<-table(popfile$t)

library(lme4)
mm0<-lmer(z~1+age+t+s+(1|ID),data=popfile)
summary(mm0)
table(popfile$C)

Camemberts<-tapply(X=popfile$C,INDEX=popfile$t,FUN=sum)#retrieve the yearly food abundance
popSize<-tapply(X=popfile$ID[which(popfile$phi==1)],INDEX=popfile$t[which(popfile$phi==1)],FUN=length)#retrieve the end of year pop sizes
tapply(X=popfile$z[which(popfile$phi==1)],INDEX=popfile$t[which(popfile$phi==1)],FUN=mean)#retrieve the end of year pop sizes
Z<-tapply(X=popfile$z,INDEX=popfile$t,FUN=mean)#retrieve the end of year pop sizes
if(length(popSize)+1==length(Camemberts)){popSize<-c(popSize,0)}

plot(Camemberts,popSize)
plot(Z,popSize)
summary(lm(popSize~1+Camemberts)
summary(lm(popSize~1+Z))
summary(lm(popSize~1+Z*Camemberts))
summary(lm(popSize~1+Z+Camemberts))


####Ellner approach####
####We must look only at adult females because only them have Annual Reproductive Success values.
females<-popfile[which(popfile$s=="F" & popfile$age>0),]

plot(density(log(females$ARS)))#very grossly gaussian distribution
#### log(ARS) ~ C+z (Hairston)

evol<-vector(length=max(females$t)-1)
ecol<-vector(length=max(females$t)-1)
for (tempus in 2:max(females$t))
{
  mARS0<-glm(ARS~1+C+z,data=females[which(females$t==(tempus-1) | females$t==tempus),],family=poisson)
  intercept<-mARS0$coefficients[1]
  dXdz<-mARS0$coefficients["z"]
  dXdk<-mARS0$coefficients["C"]
  
  Zt0<-mean(females$z[which(females$t==(tempus-1))])
  kt0<-mean(females$C[which(females$t==(tempus-1))])
  Zt1<-mean(females$z[which(females$t==tempus)])
  kt1<-mean(females$C[which(females$t==tempus)])
  
  Xtt<-mean(females$ARS[which(females$t==(tempus-1))]) # not exactly exp(intercept+Zt0*dXdz+kt0*dXdk) I do not want to account for the interaction (??) as it is problematic for XtT and XTt
  XTT<-mean(females$ARS[which(females$t==tempus)])# not exactly exp(intercept+Zt1*dXdz+kt1*dXdk) I do not want to account for the interaction (??) as it is problematic for XtT and XTt
  XtT<-exp(intercept+Zt0*dXdz+kt1*dXdk)#this is done with population mean assumption. Would be impossible on an individual basis
  XTt<-exp(intercept+Zt1*dXdz+kt0*dXdk)#this is done with population mean assumption. Would be impossible on an individual basis
  
  AX<-as.data.frame(matrix(data=c(Xtt,XTt,XtT,XTT,0,1,0,1,0,0,1,1),ncol=3))
  m1X<-(lm(V1~V2+V3,data=AX))
  
  evol[tempus-1]<-0.5*(XTT-XtT + XTt-Xtt)
  ecol[tempus-1]<-0.5*(XTT-XTt + XtT-Xtt)
}


plot(log(tapply(X=females$ARS,INDEX=females$t,FUN=mean)),ylim=c(-5,6),type="b")#log to see better what is going on at g13
points(evol/sd(evol),x=2:30-0.5,type="l",col="red")
points(ecol/sd(ecol),x=2:30-0.5,type="l",col="blue")
points(abs(evol)/abs(ecol),x=2:30-0.5,type="l",col="green")
legend(x="bottomright",legend=c("log(Noff)","Evol","Ecol","Abs(evol/ecol)"),col=c("black","red","blue","green"),lwd=2)
mean(abs(evol)/abs(ecol))


plot(ecol)
plot(abs(evol)/(abs(ecol)+abs(evol)))
mean(abs(evol)/(abs(ecol)+abs(evol)))
mean(abs(evol)/(abs(ecol)))

###############################################
########## log(ARS) ~ C+zE+zG (Ellner)#########

evol<-vector(length=max(females$t)-1)
plas<-vector(length=max(females$t)-1)
ecol<-vector(length=max(females$t)-1)
for (tempus in 2:max(females$t))
{
  mARS1<-glm(ARS~1+C+z+bvs,data=females[which(females$t==(tempus-1) | females$t==tempus),],family=poisson)
  intercept<-mARS0$coefficients[1]
  dXdz<-mARS1$coefficients["z"]
  dXdk<-mARS1$coefficients["C"]
  dXdg<-mARS1$coefficients["bvs"]
  
  Zt0<-mean(females$z[which(females$t==(tempus-1))])
  kt0<-mean(females$C[which(females$t==(tempus-1))])
  gt0<-mean(females$bvs[which(females$t==(tempus-1))])
  Zt1<-mean(females$z[which(females$t==tempus)])
  kt1<-mean(females$C[which(females$t==tempus)])
  gt1<-mean(females$bvs[which(females$t==tempus)])
  
  Xttt<-mean(females$ARS[which(females$t==(tempus-1))]) # not exactly exp(intercept+Zt0*dXdz+kt0*dXdk) I do not want to account for the interaction (??) as it is problematic for XtT and XTt
  XTTT<-mean(females$ARS[which(females$t==tempus)])# not exactly exp(intercept+Zt1*dXdz+kt1*dXdk) I do not want to account for the interaction (??) as it is problematic for XtT and XTt
  XttT<-exp(intercept+Zt0*dXdz+gt0*dXdg+kt1*dXdk)#this is done with population mean assumption. Would be impossible on an individual basis
  XtTt<-exp(intercept+Zt0*dXdz+gt1*dXdg+kt0*dXdk)#this is done with population mean assumption. Would be impossible on an individual basis
  XtTT<-exp(intercept+Zt0*dXdz+gt1*dXdg+kt1*dXdk)#this is done with population mean assumption. Would be impossible on an individual basis
  XTtt<-exp(intercept+Zt1*dXdz+gt0*dXdg+kt0*dXdk)#this is done with population mean assumption. Would be impossible on an individual basis
  XTtT<-exp(intercept+Zt1*dXdz+gt0*dXdg+kt1*dXdk)#this is done with population mean assumption. Would be impossible on an individual basis
  XTTt<-exp(intercept+Zt1*dXdz+gt1*dXdg+kt0*dXdk)#this is done with population mean assumption. Would be impossible on an individual basis

  AX<-data.frame(c(Xttt,XttT,XtTt,XtTT,XTtt,XTtT,XTTt,XTTT),c(0,0,0,0,1,1,1,1),c(0,0,1,1,0,0,1,1),c(0,1,0,1,0,1,0,1))
  names(AX)<-c("X","Z","G","K")
  m1X<-(lm(X~Z+G+K,data=AX))
  
  evol[tempus-1]<-m1X$coefficients["G"]
  plas[tempus-1]<-m1X$coefficients["Z"]
  ecol[tempus-1]<-m1X$coefficients["K"]
}
plot(log(tapply(X=females$ARS,INDEX=females$t,FUN=mean)),ylim=c(-5,6),type="b")#log to see better what is going on at g13
points(evol,x=2:30-0.5,type="l",col="red")
points(plas,x=2:30-0.5,type="l",col="purple")
points(ecol,x=2:30-0.5,type="l",col="blue")
points(abs(evol)/(abs(ecol)+abs(plas)),x=2:30-0.5,type="l",col="green")
legend(x="bottomright",legend=c("log(Noff)","Evol","Ecol","Abs(evol/ecol)"),col=c("black","red","blue","green"),lwd=2)
mean(abs(evol)/abs(ecol))
mean(abs(evol)/(abs(ecol)+abs(plas)))

plot(abs(evol)/(abs(ecol)+abs(plas)+abs(evol)),x=2:30-0.5,type="l",col="red",ylim=c(0,1))
points(abs(plas)/(abs(ecol)+abs(plas)+abs(evol)),x=2:30-0.5,type="l",col="purple")
points(abs(ecol)/(abs(ecol)+abs(plas)+abs(evol)),x=2:30-0.5,type="l",col="blue")

plot(evol)
popfile[which(popfile$t==13),]
#Conclusion: it is amazing how many fancy plots and numbers you can get out of nonsensical simulations and methods.













#quantitative genetics
library("MCMCglmm")
library("MasterBayes")
library("pedantics")
ped<-popfile[,c("ID","p1","p2")]
ped<-unique(ped[,c("ID","p1","p2")])

ped<-orderPed(ped)
names(ped)=c("id","dam","sire")
popfile$animal<-popfile$ID

#pedMass<-prunePed(ped,keep=c(levels(AllM$id[which(is.na(AllM$Weight)==F)]),levels(indinfo$id[which(indinfo$Ghost==T)])),make.base=T)

prior1 <- list(G = list(G1 = list(V = 1, nu = 0.002),G2 = list(V = 1, nu = 0.002)), R = list(V = 1,nu = 0.002))
AM0 <- MCMCglmm(z~ 1+age, random = ~animal+ID, data = popfile,pedigree=ped,prior = prior1, verbose = TRUE ,nitt = 13000, thin = 10, burnin = 3000)
summary(AM0)
plot(AM0)
autocorr(AM0$VCV)
posterior.heritability1<-AM0$VCV[,"animal"]/(AM0$VCV[,"animal"]+AM0$VCV[,"units"]+AM0$VCV[,"ID"])
plot(posterior.heritability1)
HPDinterval(posterior.heritability1,0.95)
posterior.mode(posterior.heritability1)

#just for the fun of pointing out the absurdity of our simulations:
AM1 <- MCMCglmm(z~ 1, random = ~animal+ID, data = popfile[which(popfile$age==0),],pedigree=ped,prior = prior1, verbose = TRUE ,nitt = 13000, thin = 10, burnin = 3000)
summary(AM1)
plot(AM1)
autocorr(AM1$VCV)
posterior.heritability2<-AM1$VCV[,"animal"]/(AM1$VCV[,"animal"]+AM1$VCV[,"units"]+AM1$VCV[,"ID"])
plot(posterior.heritability2)
HPDinterval(posterior.heritability2,0.95)
posterior.mode(posterior.heritability2)# of course almost 1 for new born. We should make this trait a bit more complex!
