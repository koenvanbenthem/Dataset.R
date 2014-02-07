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
#  V Survive (age based)
#  - Reproduction
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

# Growth // Sizes change, this fact is known by many - if not all
setGeneric("Grow",function(Object){standardGeneric("Grow")})

setMethod("Grow","Leprechaun",function(Object){
	Object@size<-Object@size*runif(1,1,1.2)
	return(Object)
})

# Retrieving the sex of an individual
setGeneric("Sex",function(Object){standardGeneric("Sex")})

setMethod("Sex","Leprechaun",function(Object){
	return(Object@sex)
})


############### Creating an initial population with 10 individuals
pop<-c(new("Leprechaun"))
for(i in 2:10){
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
	}
	
	#### Age+1 and growth
	for(i in ALIVE){
		pop[[i]]<-Age(pop[[i]])
		pop[[i]]<-Grow(pop[[i]])	
	}
	
	#### Reproduction  ### Not the most easy part (...)
	
	##########
	### Part dedicated to retrieving the indices of all living males and of all living females
	##########
	males<-lapply(pop,Sex)=="M" # Determine which individuals are males -- logicaly all other individuals should be females... However, this includes dead ones... Simply a list of T,T,F,F,T,F,F,....
	females<-which(!males) # Get the indices of the non-males (that is females..)
	males<-which(males)   # Get the indices of the males
	females<-intersect(females,ALIVE) # Retrieve the indices of the living(!) females
	males <-intersect(males,ALIVE) # Retrieve the indises of the living males
	
	
	##########
	# Part dedicated to breeding..  ~ But now it is weekend instead!
	##########
	
	
	### Everything should be written to a dataframe, to make sure we have all the values for ever and ever
		
}