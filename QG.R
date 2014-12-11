####### Load data and packages
require(MCMCglmm)
setwd(dir="~/GitHub/Dataset.R/Data/Simple")
#filestem <- '/home/marjob/Documents/Popdynamics/Data/'
#dat <- read.csv(paste(filestem,'TheDataSetV0.csv',sep=''),sep='\t')
dat<-read.table("TheDataSetV0.csv",header=T)
max(dat$age)

dat$p1 <- as.numeric(as.character(dat$p1))
dat$p2 <- as.numeric(as.character(dat$p2))
dat$animal <- as.factor(dat$ID)

# Create two new datasets (phenotypic data and pedigree data)
# Birth weight of individuals
exl <- dat$ID [dat$t==1 & is.na(dat$p1)]
tmp <- dat[(dat$ID %in% exl) ==F,] # exclude animals that were there in first timestep
datInd <- data.frame(
	animal = unique(tmp$animal),
	birthweight = tmp$z[duplicated(tmp$animal)==F],
	mother = tmp$p1[duplicated(tmp$animal)==F],
	father = tmp$p2[duplicated(tmp$animal)==F],
	sex = tmp$s[duplicated(tmp$animal)==F],
	C = tmp$C[duplicated(tmp$animal)==F],
	t = tmp$t[duplicated(tmp$animal)==F]
)
aa <- tapply(dat$ARS,dat$ID,sum)
tmp <- data.frame(animal=names(aa),R=as.numeric(aa))
datInd <- merge(datInd,tmp,by.x='animal')
datInd$animal <- as.factor(datInd$animal)

meantR<-tapply(datInd$R,INDEX = datInd$t,mean)
datInd$Rrelative<-0
for (time in 1:48)
{
  datInd$Rrelative[which(datInd$t==time)]<-datInd$R[which(datInd$t==time)]/meantR[time]
}


datPedi <- data.frame(
	animal=unique(dat$ID),
	mother=as.vector(tapply(dat$p1,dat$ID,function(x) x[1])),
	father=as.vector(tapply(dat$p2,dat$ID,function(x) x[1]))
)
datPedi$mother <- as.factor(datPedi$mother)
datPedi$father <- as.factor(datPedi$father)
datPedi$animal <- as.factor(datPedi$animal)

dat$age2<-(dat$age^2-mean(dat$age^2))/sd(dat$age^2) #standardize squared age

#############################
## MCMC glmm Settings #######
nitt <- 50000
thin <- 50
burnin <- 3000
varP <- var(datInd$birthweight)

# Priors for different number of random factors (all univariate)
prior1 <- list(
	G=list(
		G1=list(V=matrix(varP/2),n=1)),
	R=list(V=matrix(varP/2),n=1))

prior2<-list(
	G=list(
		G1=list(V=matrix(varP/3),n=1),
		G2=list(V=matrix(varP/3),n=1)),
	R=list(V=matrix(varP/3),n=1))

prior3<-list(
	G=list(
		G1=list(V=matrix(varP/4),n=1),
		G2=list(V=matrix(varP/4),n=1),
		G3=list(V=matrix(varP/4),n=1)),
	R=list(V=matrix(varP/4),n=1))


### Univariate
########### Birth weight
# Run animal models	
## Model 1
allModels <- list(
	model1<-MCMCglmm(birthweight~1+t,random=~animal,pedigree=datPedi,data=datInd,prior=prior1,
		nitt=nitt,thin=thin,burnin=burnin,verbose=FALSE),
	model2<-MCMCglmm(birthweight~sex,random=~animal,pedigree=datPedi,data=datInd,prior=prior1,
		nitt=nitt,thin=thin,burnin=burnin,verbose=FALSE),
	model3<-MCMCglmm(birthweight~1,random=~animal+mother,pedigree=datPedi,data=datInd,prior=prior2,
		nitt=nitt,thin=thin,burnin=burnin,verbose=FALSE),
	model4<-MCMCglmm(birthweight~sex,random=~animal+mother,pedigree=datPedi,data=datInd,prior=prior2,
		nitt=nitt,thin=thin,burnin=burnin,verbose=FALSE)
)


# Select best model
dics <- sapply(1:4,function(x) allModels[[x]]$DIC)
useModel <- allModels[[which(dics == min(dics))]]
# Choose best model and use this for analysis per timestep.

### Do it per timestep
modelsTime <- list()
for (i in unique(datInd$t)) {
	modelsTime[[i]]<-MCMCglmm(birthweight~1,random=~animal+mother,pedigree=datPedi,data=datInd[datInd$t==i,],prior=prior2,
		nitt=nitt,thin=thin,burnin=burnin,verbose=FALSE)
	print(i)
}
modelsTime <- modelsTime[-1]

# Extract parameters per timestep
allVar = sapply (modelsTime, function(x) posterior.mode(x$VCV) )#/ sum(posterior.mode(x$VCV)))

plot(unique(datInd$t),allVar[1,],type='l',lwd=2,ylim=c(0,1),xlim=c(0,50),col='red')
lines(2:50,allVar[2,],type='l',lwd=2,col='green')
lines(2:50,allVar[3,],type='l',lwd=2,col='blue')


# Explore best model
# Plot
plot(useModel$Sol) # Fixed effects
plot(useModel$VCV) # Random effects

# Calculate estimates + CI
posterior.mode(useModel$VCV) ## estimates of additive genetic, residual and other variances
HPDinterval(useModel$VCV,0.95)
mean(useModel$Sol)

# Heritabilities
heritabilities <- 
	useModel$VCV[,"animal"]/
	(useModel$VCV[,"animal"]+useModel$VCV[,"units"]) # heritability (Va / Vp)
posterior.mode(heritabilities) ## estimates of additive genetic, residual and other variances
HPDinterval(heritabilities,0.95)


###### Multivariate #######################
varP <- matrix(c(var(datInd$birthweight),0,0,var(datInd$R)),2,2)
prior <- list(
	G=list(
		G1=list(V=varP/3,n=2),
		G2=list(V=varP/3,n=2)),
	R=list(V=varP/3,n=2)
)
#priorBivariate <- list(G = list(G1 = list(V = diag(2), nu = 2),G2 = list(V = diag(2), nu = 2)), R = list(V = diag(2),nu = 2))

model<-MCMCglmm(cbind(birthweight,Rrelative)~trait-1 + trait:t,
	random=~us(trait):animal + us(trait):mother,
	rcov=~us(trait):units,
	pedigree=datPedi,data=datInd,prior=prior,
	nitt=nitt,thin=thin,burnin=burnin,verbose=TRUE,family=c('gaussian','gaussian'),pr=TRUE)
summary(model)
save(model,file='AM2.RData')
str(model)
str(model$Z@p)
# Estimate mean generation time
# Incl parents that have reproduced
#tapply(dat$age[dat$ARS>0],dat$ID[dat$ARS>0],function (x) x[1])
dat = dat[dat$phi == 1,]
allOffspring <- dat$ID[!is.na(dat$p1)]

tmp <- as.numeric()
ages <- as.numeric()
for (i in allOffspring) {
	p1 <- dat[dat$ID==allOffspring[i],]$p1[1]
	p2 <- dat[dat$ID==allOffspring[i],]$p2[1]
	t <- dat[dat$ID==allOffspring[i],]$t[1]
	tmp <- dat$age[dat$ID==p1 & dat$t==t]
	ages <- c(ages,tmp)
}

meanAgeRepr <- mean(ages)

# Use estimated genetic covA (w,z)
posterior.mode(model$VCV[,2]) * (1:50 /meanAgeRepr) + mean(dat$z[dat$t==1])

plot(1:50,tapply(dat$z,dat$t,mean),ylim=c(10,20))
lines(1:50,posterior.mode(model$VCV[,2])* (1:50 /meanAgeRepr) + mean(dat$z[dat$t==1]),col='red')

###############################################################%
###############################################################%
###########trying to retrieve the breeding values##############
###############################################################%
###############################################################%
setwd(dir="~/GitHub/Dataset.R/Data/Simple")
dat<-read.table(file="TheDataSetV0.csv",header=T)
dat$p1 <- as.numeric(as.character(dat$p1))
dat$p2 <- as.numeric(as.character(dat$p2))
dat$animal <- dat$ID

# Create two new datasets (phenotypic data and pedigree data)
# Birth weight of individuals
exl <- dat$ID [dat$t==1 & is.na(dat$p1)]
tmp <- dat[(dat$ID %in% exl) ==F,] # exclude animals that were there in first timestep
datInd <- data.frame(
  animal = unique(tmp$animal),
  birthweight = tmp$z[duplicated(tmp$animal)==F],
  mother = tmp$p1[duplicated(tmp$animal)==F],
  father = tmp$p2[duplicated(tmp$animal)==F],
  sex = tmp$s[duplicated(tmp$animal)==F],
  C = tmp$C[duplicated(tmp$animal)==F],
  t = tmp$t[duplicated(tmp$animal)==F]
)
aa <- tapply(dat$ARS,dat$ID,sum)
tmp <- data.frame(animal=names(aa),R=as.numeric(aa))
datInd <- merge(datInd,tmp,by.x='animal')
dat <- merge(dat,tmp,by.x='animal')

datInd$animal <- as.factor(datInd$animal)

meantR<-tapply(datInd$R,INDEX = datInd$t,mean)
datInd$Rrelative<-0
for (time in 1:48)
{
  datInd$Rrelative[which(datInd$t==time)]<-datInd$R[which(datInd$t==time)]/meantR[time]
}

dat$Rrelative<-dat$R/(mean(dat$R))


dat$age2<-(dat$age^2-mean(dat$age^2))/sd(dat$age^2) #standardize squared age

datPedi <- data.frame(
  animal=unique(dat$ID),
  mother=as.vector(tapply(dat$p1,dat$ID,function(x) x[1])),
  father=as.vector(tapply(dat$p2,dat$ID,function(x) x[1]))
)
datPedi$mother <- as.factor(datPedi$mother)
datPedi$father <- as.factor(datPedi$father)
datPedi$animal <- as.factor(datPedi$animal)

varP <- var(datInd$birthweight)

priorBV <- list(
  G=list(
    G1=list(V=matrix(varP/2),n=1),G2=list(V=matrix(varP/2),n=1)),
  R=list(V=matrix(varP/2),n=1))

modelBV<-MCMCglmm(birthweight~1+t,random=~animal+mother,pedigree=datPedi,data=datInd,prior=priorBV,
                 nitt=130000,thin=100,burnin=30000,verbose=TRUE,pr=TRUE)
save.image(file="AM.RData")
plot(modelBV)
autocorr(modelBV$VCV)
summary(modelBV)

mm <- data.frame(m=colMeans(modelBV$Sol),HPDinterval(modelBV$Sol))
plot(mm[order(mm$m),"m"])

BV<-modelBV$Sol[,grep(pattern = "animal*",x = colnames(modelBV$Sol))]
MV<-modelBV$Sol[,grep(pattern = "mother*",x = colnames(modelBV$Sol))]

animalID<-substr(x = colnames(BV),start = 8,stop=nchar(colnames(BV)))
motherID<-substr(x = colnames(MV),start = 8,stop=nchar(colnames(MV)))

pmBV<-data.frame(as.integer(animalID),posterior.mode(BV))
names(pmBV)<-c("ID","pBV")

pmMV<-data.frame(as.integer(motherID),posterior.mode(MV))
names(pmMV)<-c("ID","pMV")

mpmBV<-merge(x = pmBV,y = dat,by="ID",all.x=TRUE)
mpmMV<-merge(x=pmMV,y=dat,by="ID",all.x=TRUE)

head(mpmBV)
plot(mpmBV$pBV,mpmBV$bvs)
abline(a = 0,b = 1,col="red",lwd=3)
abline(lm(mpmBV$bvs ~mpmBV$pBV),col="blue",lwd=3)
par(mfrow=c(2,1))
plot(posterior.mode(BV))
plot(dat$bvs[dat$ID %in% animalID & dat$age==0])
par(mfrow=c(1,1))

summary(lm(pBV~1+t,data=mpmBV))
plot(y=mpmBV$pBV,x=mpmBV$t)
meanpBV<-tapply(mpmBV$pBV,mpmBV$t,mean)
meanrBV<-tapply(mpmBV$bvs,mpmBV$t,mean)
plot(meanpBV,type="l",ylim=c(-0.3,1))
points(meanrBV,type="l",lty=2)

DB<-meanpBV[50]-meanpBV[11]

summary(lm(pMV~1+t,data=mpmMV))
plot(y=mpmMV$pMV,x=mpmMV$t)
meanM<-tapply(mpmMV$pMV,mpmMV$t,mean)
DM<-meanM[50]-meanM[11]

meanZ<-tapply(mpmBV$z,mpmBV$t,mean)
DZ<-meanZ[50]-meanZ[11]

DR<-DZ-DM-DB

########################################################
#### Now, we should try on Z instead of Birth Weight#####
########################################################
library(lme4)

summary(lmer(z~1+t+age+age2+(1|ID),data=dat[which(dat$t>10),],REML=F))
summary(glm(z~1+t+age+age2,data=dat))

priorBV <- list(
  G=list(
    G1=list(V=matrix(varP/4),n=1),G2=list(V=matrix(varP/4),n=1),G3=list(V=matrix(varP/4),n=1)),
  R=list(V=matrix(varP/4),n=1))


modelBV<-MCMCglmm(z~1+t+age+age2,random=~animal+ID+p1,pedigree=datPedi,data=dat[which(dat$t>10),],prior=priorBV,
                  nitt=130000,thin=100,burnin=30000,verbose=TRUE,pr=TRUE)
summary(modelBV)
save(modelBV,file = "modelBV")

mm <- data.frame(m=colMeans(modelBV$Sol),HPDinterval(modelBV$Sol))
plot(mm[order(mm$m),"m"])

BV<-modelBV$Sol[,grep(pattern = "animal*",x = colnames(modelBV$Sol))]
MV<-modelBV$Sol[,grep(pattern = "p1.*",x = colnames(modelBV$Sol))]
IV<-modelBV$Sol[,grep(pattern = "ID.*",x = colnames(modelBV$Sol))]
str(modelBV$Sol)
residuals(modelBV)

animalID<-substr(x = colnames(BV),start = 8,stop=nchar(colnames(BV)))
motherID<-substr(x = colnames(MV),start = 4,stop=nchar(colnames(MV)))
IID<-substr(x = colnames(IV),start = 4,stop=nchar(colnames(IV)))

pmBV<-data.frame(as.integer(animalID),posterior.mode(BV))
names(pmBV)<-c("ID","pBV")

pmMV<-data.frame(as.integer(motherID),posterior.mode(MV))
names(pmMV)<-c("ID","pMV")

pmIV<-data.frame(as.integer(IID),posterior.mode(IV))
names(pmIV)<-c("ID","pIV")

mpmBV<-merge(x = pmBV,y = dat,by="ID",all.x=TRUE)
mpmMV<-merge(x=pmMV,y=dat,by="ID",all.x=TRUE)
mpmIV<-merge(x=pmIV,y=dat,by="ID",all.x=TRUE)

head(mpmBV)
plot(mpmBV$pBV,mpmBV$bvs)
abline(a = 0,b = 1,col="red",lwd=3)
abline(lm(mpmBV$bvs ~mpmBV$pBV),col="blue",lwd=3)
par(mfrow=c(2,1))
plot(posterior.mode(BV))
plot(dat$bvs[dat$ID %in% animalID & dat$age==0])
par(mfrow=c(1,1))

summary(lm(pBV~1+t,data=mpmBV))
plot(y=mpmBV$pBV,x=mpmBV$t)
meanpBV<-tapply(mpmBV$pBV,mpmBV$t,mean)
meanrBV<-tapply(mpmBV$bvs,mpmBV$t,mean)
plot(meanpBV,type="l",ylim=c(-0.3,1))
points(meanrBV,type="l",lty=2)

DB<-meanpBV[50]-meanpBV[11]


summary(lm(pMV~1+t,data=mpmMV))
plot(y=mpmMV$pMV,x=mpmMV$t)
meanM<-tapply(mpmMV$pMV,mpmMV$t,mean)
DM<-meanM[50]-meanM[11]

meanI<-tapply(mpmIV$pIV,mpmIV$t,mean)
DI<-meanI[50]-meanI[11]

meanZ<-tapply(mpmBV$z,mpmBV$t,mean)
DZ<-meanZ[50]-meanZ[11]

DR<-DZ-DM-DB-DI
DR
DM 
DB 
DI

#Cleaner: now accounting for uncertainty
meanT<-tapply(mpmBV$t,mpmBV$ID,mean)
lmBV<-as.mcmc(apply(BV,MARGIN = 1,function(x){coef(lm(x~1+meanT))[2]}))
plot(lmBV)
posterior.mode(lmBV*39);HPDinterval(lmBV*39)# estimated slope and CI

meanTID<-tapply(mpmBV$t,mpmBV$ID,mean)

meanTID<-meanTID[paste("ID.",names(meanTID),sep="")%in%colnames(IV)]
mean(paste("ID.",names(meanTID),sep="")==colnames(IV))
lmIV<-as.mcmc(apply(IV,MARGIN = 1,function(x){coef(lm(x~1+meanTID))[2]}))
plot(lmIV)
posterior.mode(lmIV*39);HPDinterval(lmIV*39)# estimated slope and CI

meanTmother<-tapply(mpmBV$t,mpmBV$p1,mean)
meanTmother<-meanTmother[paste("p1.",names(meanTmother),sep="")%in%colnames(MV)]
mean(paste("p1.",names(meanTmother),sep="")==colnames(MV))
lmMV<-as.mcmc(apply(MV,MARGIN = 1,function(x){coef(lm(x~1+meanTmother))[2]}))
plot(lmMV)
posterior.mode(lmMV*39);HPDinterval(lmMV*39)# estimated slope and CI

lmRV<-DZ-lmMV*39-lmIV*39-lmBV*39
plot(lmRV)
posterior.mode(lmRV)

posteriorDistriOfChange<-list(lmBV*39,lmIV*39,lmMV*39,lmRV)
names(posteriorDistriOfChange)<-c("BreedingValues","Identity","Mother","Residual")
posteriorModeOfChange<-list(posterior.mode(lmBV)*39,posterior.mode(lmIV)*39,posterior.mode(lmMV)*39,posterior.mode(lmRV))
names(posteriorModeOfChange)<-c("BreedingValues","Identity","Mother","Residual")
posteriorHPDOfChange<-list(HPDinterval(lmBV)*39,HPDinterval(lmIV)*39,HPDinterval(lmMV)*39,HPDinterval(lmRV))
names(posteriorHPDOfChange)<-c("BreedingValues","Identity","Mother","Residual")

ResultsQG<-list(posteriorModeOfChange,posteriorHPDOfChange,posteriorDistriOfChange)
names(ResultsQG)<-c("mode","CI","CompleteDistri")
save(ResultsQG,file = "ResultsQG.RData")

ResultsQG$mode
save.image(file = "AM.RData")

##### Bivariate#####
library(lme4)
summary(lmer(Rrelative~1+age+age2+(1|t)+(1|ID),data=dat))
summary(lm(Rrelative~1+t+age+age2,data=dat))

b <- 1e+10
B1_V <- diag(5)*b; B1_mu <- rep(0,5)
B1_V[2,2]<-0.0001
B1_mu[2]<-1

priorMV <- list(
  G=list(
    G1=list(V=diag(2),n=2),
    G2=list(V=diag(2),n=2),
    G3=list(V=diag(1),n=1)),
  R=list(V=matrix(data=c(1,0,0,0.00001),nrow = 2),n=2,fix=2),
  B=list(V=B1_V, mu=B1_mu)
)
#priorBivariate <- list(G = list(G1 = list(V = diag(2), nu = 2),G2 = list(V = diag(2), nu = 2)), R = list(V = diag(2),nu = 2))

modelMV<-MCMCglmm(cbind(z,Rrelative)~trait-1 + at.level(trait,c(1)):(age+age2+t),
                random=~us(trait):animal + us(trait):ID + us(at.level(trait,c(1))):p1,
                rcov=~idh(trait):units,
                pedigree=datPedi,data=dat[which(dat$t>10),],prior=priorMV,
                nitt=150000,thin=100,burnin=50000,verbose=TRUE,family=c('gaussian','gaussian'),pr=TRUE)

save(modelMV,file='AMMV.RData')
summary(modelMV)

meanrBV<-tapply(dat$bvs,dat$t,mean)
plot( meanrBV,xlim=c(-10,60),ylim=c(-2,2))
abline(a = meanrBV[10]-0.02297*10,b=0.02297,lwd=3)
abline(a = meanrBV[10]-0.02297*10,b=-0.1136)
abline(a = meanrBV[10]-0.02297*10,b=0.1772)


priorFitness <- list(  G=list(
  G1=list(V=0.1,n=10),
  G2=list(V=7,n=0.2)),
  R=list(V=0.00001,n=0.002,fix=1),
  B=list(mu=1,V=0.0001))
Starting<-list(G=list(0.1,4.5),R=list(0.00001))

modelFitness<-MCMCglmm(Rrelative~1,random=~animal+ID,pedigree=datPedi,
                  data=dat[which(dat$t>10),],prior=priorFitness,start=Starting,
                  nitt=13000,thin=10,burnin=0,verbose=TRUE,pr=TRUE)
summary(modelFitness)
plot(modelFitness)

var(dat$R)
summary(lmer(Rrelative~1+age+(1|ID),data=dat))
summary(lmer(R~1+(1|ID),data=dat))
