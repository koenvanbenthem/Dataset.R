
#######################################################
########### Analyze simulated data on Leprechauns #####
#################### June 2014 ########################
################# Marjolein Bruijning #################
#######################################################


filestem='/home/marjob/Documents/Popdynamics/'

# Load all data files & description
folder <- 'Test4' # folder containing csvs
filenames <- list.files(path=paste(filestem,"Data/",folder,sep=""))
filenames <- sort(filenames)
allData <- lapply(paste(filestem,"Data/",folder,"/",filenames,sep=""), read.csv,sep="")
descrip <- read.csv(paste(filestem,"Data/conv.csv",sep=""),sep=",")
descrip$filename <- substr(descrip$filename,12,100) # remove first 11 characters
descrip <- descrip[order(descrip$filename),] # same ordre as filenames

# Remove populations with <n individuals
n <- 50
exl <- c(which(unlist(lapply(allData,function(x)dim(x)[1])) < n))
if (length(exl) > 0) {
	allData <- allData [-exl]
	descrip <- descrip[-exl,]
}

############################################
################# Run! #####################
############################################

# Analyze all datasets
allResults <- lapply(allData, analyzeAll)

### Select datasets to analyze
## Environmental variation
inc <- which(descrip[,2] == 0.1 & descrip[,3] == 0 & descrip[,4] == 0 & descrip[,5] == 0 & descrip[,6] == 0 )
allResults <- lapply(allData[inc], analyzeAll)
plot(descrip$sdc[inc],unlist(lapply(allResults,function(x) x$camembert$sd)), main='Environmental variation',
	type='b',xlab='Model input',ylab='Estimated value from simulated data')

## Viability selection
inc <- which(descrip[,3] == 0 & descrip[,4] == 0 & descrip[,5] == 0 & descrip[,6] == 0 & descrip[,7] == 1000 )
allResults <- lapply(allData[inc], analyzeAll)
plot(descrip$SS1[inc],unlist(lapply(allResults,function(x) x$survSel[1,2])),main='Viability selection',
	type='p',xlab='Model input',ylab='Estimated value from simulated data')

## Fertility selection (ARS)
inc <- which(descrip[,2] == 0.1 & descrip[,3] == 0 & descrip[,5] == 0 & descrip[,6] == 0 & descrip[,7] == 1000 )
allResults <- lapply(allData[inc], analyzeAll)
plot(descrip$fs2[inc],unlist(lapply(allResults,function(x) x$reproSel[1,2])),main='Fertility selection',
	type='b',xlab='Model input',ylab='Estimated value from simulated data')



###########################################
################ Functions ################
###########################################

analyzeAll = function (dat) {
	# Prepare dataset
	dat$p1 <- as.numeric(as.character(dat$p1))
	dat$p2 <- as.numeric(as.character(dat$p2))
	dat$animal <- dat$ID
	results <- list ()
	
	## Environmental variation
	results$camembert <- data.frame(sd=sd(tapply(dat$C,dat$t,sum)), # sd total Camembert abundance
		mean=mean(tapply(dat$C,dat$t,sum))) # mean total Camembert abundance

	## Differences between mother and offspring in birth weight
	allID=unique(dat$ID)
	allID = allID[(allID %in% dat$ID[dat$t==1]) ==F] # Exclude first cohort individuals
	
	if (length(allID) > 1 ) {
		z=data.frame(offspring=NA,mother=NA,father=NA,allID=allID)

		for (i in 1:length(allID)) {
			inc=dat$ID == allID[i]
			z$offspring[i] <- dat$z[inc][1]
			
			z$mother[i] <- dat$z[dat$ID == dat$p1[inc][1]][1]
			z$father[i] <- dat$z[dat$ID == dat$p2[inc][1]][1]
		}

		z$meanParent <- (z$mother + z$father)/2
		if (dim(z)[1] > 2) {
			results$birthweight = data.frame(
				h2=summary(lm(z$offspring~z$meanParent))$coefficients[2,1],
				Dbar=mean(z$offspring - z$meanParent,na.rm=T),
				DbarSd=sd(z$offspring - z$meanParent,na.rm=T)
			)
		} else {results$birthweight <- NA}
	} else {
		results$birthweight <- NA
	}
	
	#### (Plasticity in) Growth
	allID <- unique(dat$ID)
	allID <- allID[(allID %in% dat$ID[dat$t==1]) ==F] # Exclude first cohort individuals

	if (length(allID) > 1 ) {
		#plot(110,110,xlim=c(0,20),ylim=c(0,20))
		growth <- data.frame(growth=NA,z=NA,IDs=NA)
		#plot(110,110,xlim=c(0,20),ylim=c(0,20))
		for (i in 1:length(allID)) {
			inc <- dat$ID == allID[i] & dat$phi==1 # Exlude dead individuals
			#lines(dat$t[inc],dat$z[inc])
			tmp <- data.frame (
				growth=dat$z[inc][-1] / dat$z[inc][-length(dat$z[inc])], # Proportional growth
				z=dat$z[inc][-length(dat$z[inc])],
				IDs=dat$ID[inc][-1]
			)
			growth <- rbind(growth,tmp)
		}
		growth <- growth[-1,]

		#hist(growth$growth)
		if (dim(growth)[1] > 0) {
			results$withtinPlasticity <- data.frame(
				mingrowth=min(growth$growth,na.rm=T),
				maxgrowth=max(growth$growth,na.rm=T)
			)
		} else {results$withinPlasticity <- NA}
	} else {
		results$withinPlasticity <- NA
	}
	##### Selection
	# Exluce all dead individuals from now on
	dat=dat[dat$phi == 1,]
	dat$surv = 1

	for (i in 1:length(unique(dat$ID))) {
		inc = dat$ID == unique(dat$ID)[i]
		dat$surv[inc] [length(dat$surv[inc])] = 0 # Replace last alive census by 0
	}

	#par(mfrow=c(3,3))
	dat$sizeDev <- dat$z-(10*1.045^dat$age)
	## Viability selection
	ages=c(0:2)
	results$survSel = data.frame(age=ages,coeff=NA)
	for (i in 1:length(ages)) {
		inc = dat$age == ages[i]
		plot(dat$z[inc],dat$surv[inc],main=paste('Age: ',unique(dat$age)[i],sep=''))
		mod=glm(surv~sizeDev + sizeDev^2,data=dat[inc,],family='binomial')
		results$survSel$coeff[i] = mod$coefficients[2]
		lines(1:15,predict(mod,newdata=data.frame(z=1:15),type='response'),col='red')
	}

	## Fertility selection
	#par(mfrow=c(3,3))
	ages=1:3
	results$reproSel = data.frame(age=ages,coeff=NA)
	for (i in 1:length(ages)) {
		inc = dat$age == ages[i]
		#plot(dat$z[inc],dat$ARS[inc],main=paste('Age: ',unique(dat$age)[i],sep=''))
		mod=glm(ARS~z,data=dat[inc,],family='poisson')
		results$reproSel$coeff[i] = mod$coefficients[2]
		#lines(1:15,predict(mod,newdata=data.frame(z=1:15),type='response'),col='red')
	}
	print(Sys.time())
	return(results)
}

