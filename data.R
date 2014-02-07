cat("\014")

setwd(dir="C:/Users/Thimothee Admin/Documents/thesis/Transitions Life Histories/CppCodeSimul/")
library(lme4)

NormalROracle<-as.data.frame(matrix(data=NA,nrow=0,ncol=61))
names(NormalROracle)<-c("SamplingSeed","PoissonSimulNumber","NeutralSimulNumber","AdultSurvival","JuvenileSurvival","IndVarSurv","IndVarRepro"
                        ,"CorReproSurv","MeanRepro","AgeMax","StudyLength","IndInCohort"
                        ,"SSLRTreject","RRLRTreject","ViS2.5%","ViS50%","ViS97.5%","ViR2.5%","ViR50%","ViR97.5%"
                        ,"SurvReproreject","SurvAgereject","ReproAgereject","PluriReproreject"
                        ,"SR2.5%","SR50%","SR97.5%","SA2.5%","SA50%","SA97.5%","RA2.5%","RA50%","RA97.5%","PR2.5%","PR50%","PR97.5%"
                        ,"RQ2.5%","RQ50%","RQ97.5%","SQ2.5%","SQ50%","SQ97.5%","ReproQreject","SurvQreject"
                        ,"ReproQR2_2.5%","ReproQR2_50%","ReproQR2_97.5%","SurvQR2_2.5%","SurvQR2_50%","SurvQR2_97.5%"
                        ,"SRQR2.5%","SRQR50%","SRQR97.5%","SRQS2.5%","SRQS50%","SRQS97.5%","SRQRreject","SRQSreject"
                        ,"SRQR2_2.5%","SRQR2_50%","SRQR2_97.5%")
##### In this file we will try to come up with some things ####
###############################################################

########################
####### Type of program
## Written in R
## Object Oriented
#######################

#######################
# Main object "Leprechaun": an individual of our species of interest, the Irish Leprechaun, small in size, but with great powers. To make it more French for Timothee, we assume that their favourite food is camembert.
#
# The object contains the following values

# Static [these numbers do not change after initialisation]
#  V ID (unique identifier number)																	[integer]
#  V pID (two numbers referring to the parents of the Leprechaun, if none: NA)						[vector of two integers]
#  V Year of birth (timestep at which the individual was born)										[integer]
#  - Genome (?) (two vectors of length N coding for both chromosomes of N loci in the genome.)		[two vectors of N integers]

# Dynamic [these numbers do change after initialisation]
#  V alive (boolean, true/false)																	[boolean]
#  V age (timesteps since birth)																	[integer]
#  V size																							[double]

# The object contains the following functions (for now)
#  - Grow
#  - Survive
#######################

####################################################
################ GLOBAL VARIABLES AND COUNTERS #####
####################################################

################ Counter for the IDs ###############
CID<-as.integer(1) 


################ Counter for the current year ######
YR<-0 

####################################################
############### Definition of the class ############
####################################################

setClass(
	Class="Leprechaun",
	representation=representation(
		ID = "integer",
		pID = "integer",
		age = "integer",
		Birth = "integer",
		alive = "logical",
		size = "numeric",
		sex = "character"
	)

)

############### Definition of the basic methods (for printing to the screen and initialisation)
###############################################################################################
setMethod("show","Leprechaun",
	function(object){
		cat(object@ID,"\t",object@size,"\t",object@age,"\t",object@sex,"\t(",object@pID[1],",",object@pID[2],")\t",object@Birth,"\t",object@alive,"\n",sep="")
	}
)

setMethod("initialize","Leprechaun",function(.Object,parent1,parent2){
	if(missing(parent1)){parent1<-NA; weight1<-5+2*runif(1)}else{weight1<-pop[[parent1]]@size}
	if(missing(parent2)){parent2<-NA; weight2<-5+2*runif(1)}else{weight2<-pop[[parent2]]@size}
	.Object@age<-as.integer(0)
	.Object@ID<-CID
	.Object@pID<-c(as.integer(parent1),as.integer(parent2))
	.Object@Birth<-as.integer(YR)
	.Object@alive<-TRUE
	.Object@size<-0.5*weight1+0.5*weight2

	if(runif(1)>0.5){.Object@sex<-'F'}else{.Object@sex<-'M'}

	CID<<-as.integer(CID+1)
	
	return(.Object)
})

################### Definition of more biologically relevant methods (e.g. survival)
####################################################################################

# Implementing the famous bathtub, ages 1 to 20
bathtub<-function(age){
	p<-0.6*exp(-age/4)+(-1+exp(age*log(2)/20))
	p[p>1]<-1
	return(p)
}

# Applying the bathtub in a surival function
setGeneric("Surv",function(Object){standardGeneric("Surv")})

setMethod("Surv","Leprechaun",function(Object){
	
	if(runif(1)>bathtub(Object@age)){
		Object@alive<-FALSE
		ALIVE<<-ALIVE[ALIVE!=Object@ID]
	}
	
	return(Object)
	
})

# Simple function, simply adds 1 to the age
setGeneric("Age",function(Object){standardGeneric("Age")})

setMethod("Age","Leprechaun",function(Object){
	Object@age<-as.integer(Object@age+1)
	return(Object)
})


############### Creating an initial population with 10 individuals
pop<-c(new("Leprechaun"))
for(i in 2:100){
	pop<-c(pop,new("Leprechaun"))
}

############### List of living individuals [their indices], this will save time later, because dead individuals are not looped over
ALIVE<-1:length(pop)

############### The start of time
for(YR in 1:10){
	cat("\nAt the beginning of year:",YR,"\nThere are:",length(ALIVE),"Leprechauns\n-----------------\n")
	
	#### Survival
	for(i in ALIVE){
		pop[[i]]<-Surv(pop[[i]])
		pop[[i]]<-Age(pop[[i]])
		cat(i,"\t")
	}
	
	
	#### Growth
	
	
	#### Reproduction
		
}