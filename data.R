##### In this file we will try to come up with some things ####
###############################################################

########################
####### Type of program
## Written in R
## Object Oriented
#######################

#######################
# Main object "Lepricorn": an individual of our species of interest, the Irish lepricorn, small in size, but with great powers. To make it more French for Timothee, we assume that their favourite food is camembert.
#
# The object contains the following values

# Static [these numbers do not change after initialisation]
#  V ID (unique identifier number)																	[integer]
#  V pID (two numbers referring to the parents of the lepricorn, if none: NA)						[vector of two integers]
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
	Class="Lepricorn",
	representation=representation(
		ID = "integer",
		pID = "integer",
		age = "integer",
		Birth = "integer",
		alive = "logical",
		size = "numeric"
	)

)

############### Definition of the basic methods (for printing to the screen and initialisation)
###############################################################################################
setMethod("show","Lepricorn",
	function(object){
		cat(object@ID,"\t",object@size,"\t",object@age,"\t(",object@pID[1],",",object@pID[2],")\t",object@Birth,"\t",object@alive,"\n",sep="")
	}
)

setMethod("initialize","Lepricorn",function(.Object,parent1,parent2){
	if(missing(parent1)){parent1<-NA; weight1<-5+2*runif(1)}else{weight1<-pop[[parent1]]@size}
	if(missing(parent2)){parent2<-NA; weight2<-5+2*runif(1)}else{weight2<-pop[[parent2]]@size}
	.Object@age<-as.integer(0)
	.Object@ID<-CID
	.Object@pID<-c(as.integer(parent1),as.integer(parent2))
	.Object@Birth<-as.integer(YR)
	.Object@alive<-TRUE
	.Object@size<-0.5*weight1+0.5*weight2
	CID<<-as.integer(CID+1)
	
	return(.Object)
})

############### Creating an initial population with 10 individuals
pop<-c(new("Lepricorn"))
for(i in 2:10){
	pop<-c(pop,new("Lepricorn"))
}