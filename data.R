##### In this file we will try to come up with some things ####
###############################################################

########################
####### Type of program
## Written in R
## Object Oriented
#######################

#######################
# Main object "Leprechaun": an individual of our species of interest, the Irish Leprechaun, small in size, but with great powers. To make it more French for Timothee, we assume that their favourite food is camembert. Leprechaun is not very choosy, and mates completely random.
##########################################

# The object contains the following values

# Static [these numbers do not change after initialisation]
#  V ID (unique identifier number) [integer]
#  V pID (two numbers referring to the parents of the Leprechaun, if none: NA) [vector of two integers]
#  V Year of birth (timestep at which the individual was born)	[integer]
#  (- Genome (?) (two vectors of length N coding for both chromosomes of N loci in the genome.) [two vectors of N integers]
#  - Heritable phenotypic trait value of interest (z) (e.g. birth weight) (changes/constant through life depending on trait)
#  - (Possibly: breeding value A))
#  - Rather simulate physically independent loci, otherwise we need to simulate recombination on the chromosomes. Over short time periods, 100 recombinations fragments (ie independent loci) sounds realistic.
#  - Explicit coding of traits by many independent diploid loci. The simplest model: z = mean + sum_loci(a1_locus + a2_locus) + environment. One can add explicit dominance and epistasis, as well as interactions with environment. 
#  - A possibility is to draw the a of the different alleles from a N(0,V). Each locus can have a different V, and thus a different importance. 
#  - A large number of loci (>20) will give easily patterns expected from quantitative genetics. We can draw randomly the number of loci per trait. 
#  - My main concern at the moment is the initialisation of the genetic diversity: it will be hard to avoid a fast decline of diversity at the beginning. One possibility is using neutral expectations of diversity (n-coalescent or Ewens distribution)
#  - We probably do not need mutations if we consider a population over no more than some tens of generations.
#  V Sex (M/F). Could be genetically determined by on locus, thus allowing random fluctuations of sex ratio and thus population structure.

# Dynamic [these numbers do change after initialisation]
#  V alive (boolean, true/false) [boolean]
#  V age (a) (timesteps since birth) [integer]
#  - stage (possibly instead of/in addition to age) (juvenile, adult, etc)
#  V size (x)

####################
###### Relations that need to be defined between size (x), age (a), trait (z) and vital rates
# (- Possibly include population density d in functions)
# V Survival(a,x,z,d) (logistic function) 
# - Growth (a,x,z,d) (either transition probability to next stage or absolute growth)
# - Reproduction probability: p_repr(a,x,z,d) (logistic function)
# - Number of offspring: n_offspring(a,x,z,d) (poisson distribution)
# - Offspring size x distribution: x_offspring(x_mother,z_mother,x_father,z_father) (not required if working with stages) (prability density function)
# - Offspring trait z distribution: z_offspring(x_mother,z_mother,x_father,z_father,A_mother,A_father) + rnorm(0,V_E) (prability density function)
# - For survival and reproduction, we can consider functions of the form W=exp(-((z-Optimum)^2)/(2*width)), that is, gaussian stabilizing selection. We could then let the Optimum shift. 

####################
####### Environmental aspects
# - Extra variation in z due to unexplained factors (i.e. V_E)) (assumed constant over time and constant accross all individuals)
# - Changes in the environment affecting survival(x,z), growth(x,z), p_repr(x,z), n_offspring(x,z), f_offspring(x_mother,z_mother) (i.e. changing selection)

########################
########## 'Settings' of ancestral population from which simulation can start:
# - n start individuals
# - Start trait z distribution
# - Start age (/stage) distribution
# - Assign sexes to individuals
# - (Possibly: start a values)

########## Perform in each time step the following actions over all alive individuals at t=0
## In following order:
# survival(x,z)
### For those who survive:
# growth(x,z)
# p_repr(x,z)
### For those who reproduce:
# Random mating between all reproductive males and females
# n_offspring(x,z)
# x_offspring(x_mother,z_mother,x_father,z_father)
# z_offspring(x_mother,z_mother,x_father,z_father,a_mother,a_father,m) + rnorm(0,V_E)
# Random sex assigned to offspring
# Offspring added to population
### End of timestep

####################################################
################ GLOBAL VARIABLES AND COUNTERS #####
####################################################

################ Random seed #######################
set.seed(12)

################ Counter for the IDs ###############
CID<-as.integer(1) 


################ Counter for the current year ######
YR<-0 

################ Base population parameters ########

MeanBirthSize<-10
lowBoundGrowth<-1 # minimal growth rate
highBoundGrowth<-1.1 # maximal growth rate
meanGrowth<-mean(runif(10000,lowBoundGrowth,highBoundGrowth)) # needed to make survival ~ size relative on age
MeanRepro<-2

################ Environmental parameters ##########
MeanCamembert<-5000
SDCamembert<-1000


############### Genetic determinism ################
dominance<-1 # for additive effects only, must be 0
overdominance<-0 # non-null values generate overdominance

nbLoci<-10 #number of loci controling the trait phenotype
nbAlleles<-10 #number of existing alleles per loci

gvalues<-array(data=NA,dim=c(nbAlleles,nbAlleles,nbLoci),dimnames=list(paste("A",1: nbAlleles,sep=""),paste("A",1: nbAlleles,sep=""),paste("L",1:nbLoci,sep=""))) # Initialising a matrix that will contain the genotypic effects on the/a trait

for(L in 1:nbLoci)
{
  # Setting the effects for the homozygotes [all loci]
  effect<-abs(rnorm(n=1,mean=0,sd=1))# alter the locus importance in a realistic way (many small-effect loci, few major loci)
  diag(gvalues[,,L])<-2*rnorm(n=dim(gvalues)[1],mean=0,sd=effect)
  
  
  # Setting the effects for the heterozygotes
  for(A in 1:(nbAlleles-1))# loop for off-diagonal = heterozygotes (additive and dominance effects)
  {
    for (D in (A+1):nbAlleles)
    {
      d<-dominance*runif(n=1,min=-0.5-overdominance,max=0.5+overdominance)
      gvalues[A,D,L]<-(0.5-d)*gvalues[A,A,L]+(0.5+d)*gvalues[D,D,L] # mean of additive effects + dominance, over diagonal
      gvalues[D,A,L]<-(0.5-d)*gvalues[A,A,L]+(0.5+d)*gvalues[D,D,L] # the same below diagonal    
    }
  }
}
############### Selection parameters ###############
SurvivalSelection1<-0.1 #linear coefficient on a logit scale for Survival ~ ... + size +size^2
SurvivalSelection2<-(-0.01) #quadratic coefficient on a logit scale for Survival ~ ... + size + size^2; negative value=balancing selection

fertilitySelection1<-0.1 #linear coefficient on a log scale for reproduction ~ ... + size + size^2
fertilitySelection2<-(-0.01) #quadratic coefficient on a log scale for reproduction ~ ... + size + size^2; negative value=balancing selection

camembertSelection<-0.1
survivalPenaltyForRepro<-0
####################################################


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
    camemberts = "integer",
		sex = "character",
		DNA = "matrix",
    bvs = "numeric",
    ARS = "integer"
	)

)

############### Definition of the basic methods (for printing to the screen and initialisation)
###############################################################################################
setMethod("show","Leprechaun",
	function(object){
		cat(object@ID,"\t",object@size,"\t",object@camemberts,"\t",object@age,"\t",object@sex,"\t(",object@pID[1],",",object@pID[2],")\t",object@Birth,"\t",object@alive,"\n",sep="")
	}
)

setMethod("initialize","Leprechaun",function(.Object,parent1,parent2){
	.Object@DNA<-matrix(NA,nrow=2,ncol=nbLoci)
	if(missing(parent1)){
		parent1<-NA
		#weight1<-MeanBirthSize+2*runif(1)
		.Object@DNA[1,]<-floor(runif(nbLoci,min=1,max=nbAlleles+1))
	}else{
			weight1<-pop[[parent1]]@size
			.Object@DNA[1,]<-pop[[parent1]]@DNA[cbind(floor(runif(n=nbLoci,min=1,max=3)),1:nbLoci)]
	}
	if(missing(parent2)){
		parent2<-NA
		#weight2<-MeanBirthSize+2*runif(1)
		.Object@DNA[2,]<-floor(runif(nbLoci,min=1,max=nbAlleles+1))
	}else{
		#weight2<-pop[[parent2]]@size
		.Object@DNA[2,]<-pop[[parent2]]@DNA[cbind(floor(runif(n=nbLoci,min=1,max=3)),1:nbLoci)]
	}
	.Object@age<-as.integer(0)
	.Object@ID<-CID
	.Object@pID<-c(as.integer(parent1),as.integer(parent2))
	.Object@Birth<-as.integer(YR)
	.Object@alive<-TRUE
	#.Object@size<-0.5*weight1+0.5*weight2
  BreedingValueSize<-0
  for (Locus in 1:nbLoci)#take the mean of genetic values
    {
      BreedingValueSize<-BreedingValueSize+(gvalues[ .Object@DNA[1,Locus], .Object@DNA[2,Locus], Locus]/nbLoci)
    }
  .Object@bvs<-BreedingValueSize
	size<-MeanBirthSize+BreedingValueSize
  .Object@size<-rnorm(n=1,mean=size,sd=1)
  .Object@camemberts<-as.integer(0)
  
  .Object@ARS<-as.integer(0)#annual reproductive success
  
	if(runif(1)>0.5){.Object@sex<-'F'}else{.Object@sex<-'M'}

	CID<<-as.integer(CID+1)
	
	return(.Object)
})


################### Definition of more biologically relevant methods (e.g. survival)
####################################################################################

# Implementing the famous bathtub, ages 1 to 20
bathtub<-function(age,size){
  p<-0.6*exp(-age/4)+(-1+exp(age*log(2)/20))
  p[p>1]<-1
  return(p)
}

# Incorporate the effect of size to the survival function
sizeSurvival<-function(age,size,camemberts){
  sizedeviation<-size-(MeanBirthSize*meanGrowth^age)
  p<-bathtub(age)
  if(p<1)#because size does not prevent animals of maximal age to die out
    {
      plogit<-log(p/(1-p))
      Philogit<-plogit-SurvivalSelection1*sizedeviation-SurvivalSelection2*sizedeviation^2
      p<-exp(Philogit)/(1+exp(Philogit))
      if(camemberts<100)
        {
          p<-p*(camemberts/100)
        }
    }
  return(p)
}

# Applying the bathtub in a surival function
setGeneric("Surv",function(Object){standardGeneric("Surv")})

setMethod("Surv","Leprechaun",function(Object){
	
	if(runif(1)>sizeSurvival(Object@age,Object@size,Object@camemberts)){
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
	Object@size<-Object@size*runif(1,lowBoundGrowth,highBoundGrowth)
	return(Object)
})

# Retrieving the ID of an individual
setGeneric("IDretrieve",function(Object){standardGeneric("IDretrieve")})

setMethod("IDretrieve","Leprechaun",function(Object){
  return(Object@ID)
})

# Retrieving the size of an individual
setGeneric("Size",function(Object){standardGeneric("Size")})

setMethod("Size","Leprechaun",function(Object){
  return(Object@size)
})

# Retrieving the sex of an individual
setGeneric("Sex",function(Object){standardGeneric("Sex")})

setMethod("Sex","Leprechaun",function(Object){
  return(Object@sex)
})

# Calculating the number of offspring for a females
setGeneric("Num_off",function(Object){standardGeneric("Num_off")})

setMethod("Num_off","Leprechaun",function(Object){
  sizeDeviation<-Object@size-MeanBirthSize*meanGrowth^Object@age    
  lambda<-exp(log(MeanRepro)+fertilitySelection1*sizeDeviation+fertilitySelection2*sizeDeviation^2+camembertSelection*(sqrt(Object@camemberts)-survivalPenaltyForRepro))
  repro<-rpois(n=1,lambda=lambda)
  Object@ARS<-as.integer(repro)
	return(Object)
})

# Function for camembert attributions 
setGeneric("Food",function(Object){standardGeneric("Food")})

setMethod("Food","Leprechaun",function(Object){
  Object@camemberts<-as.integer(podium[i])
  return(Object)
})


############### Creating an initial population with 10 individuals
pop<-c(new("Leprechaun"))
for(i in 2:10){
	pop<-c(pop,new("Leprechaun"))
}

############### List of living individuals [their indices], this will save time later, because dead individuals are not looped over
ALIVE<-1:length(pop)

filename<-"pop.csv"
cat("t\tID\tz\tbvs\tC\ts\tARS\tage\tp1\tp2\tphi",file=filename,append=FALSE)

############### The start of time
for(YR in 1:30){
  camembert<-abs(round(rnorm(n=1,mean=MeanCamembert,sd=SDCamembert),digits=0)) # ressources for year YR
	cat("\nAt the beginning of year:",YR,"\nThere are:",length(ALIVE),"Leprechauns\n-----------------\n")
	cat("ALIVE:",ALIVE,"\n")
	cat("Current ID:",CID,"\n")
  cat("\n-----------------\nThe camembert production is:",camembert,"\n")
	#### Competition for ressources
  SizesAlive<-as.numeric(lapply(pop,Size))[ALIVE]
  #HunterQualities<-rep(x=1/length(ALIVE),length(ALIVE))# here completely random. prob can introduce quality for competition
  HunterQualities<-exp(SizesAlive)/sum(exp(SizesAlive))# here larger individual have a strong advantage in the competition
  podium<-table(sample(as.character(ALIVE),size=camembert,replace=T,prob=HunterQualities)) 
  for(i in 1:length(podium)){
    pop[[as.integer(names(podium))[i]]]<-Food(pop[[as.integer(names(podium))[i]]])
  }
  
  #### Survival
  FATUM<-ALIVE
	for(i in ALIVE){
		pop[[i]]<-Surv(pop[[i]])
	}
	DEAD<-setdiff(FATUM,ALIVE)
	#### Age+1 and growth
	for(i in ALIVE){
		pop[[i]]<-Age(pop[[i]])
		pop[[i]]<-Grow(pop[[i]])	
	}
	
	#### Reproduction  ### Not the easiest part (but necessary if we want to allow the population to grow)
	
	##########
	### Part dedicated to retrieving the indices of all living males and of all living females
	##########
	males<-lapply(pop,Sex)=="M" # Determine which individuals are males -- logicaly all other individuals should be females... However, this includes dead ones... Simply a list of T,T,F,F,T,F,F,....
	females<-which(!males) # Get the indices of the non-males (that is females..)
	males<-which(males)   # Get the indices of the males
	females<-intersect(females,ALIVE) # Retrieve the indices of the living(!) females
	males <-intersect(males,ALIVE) # Retrieve the indises of the living males
	#cat(females)
	
	##########
	# Part dedicated to breeding.. 
	##########
	
	# We take a female based approach: we determine for each females 
	from<-CID
	for(i in females){
		pop[[i]]<-Num_off(pop[[i]])
    Noffs<-pop[[i]]@ARS
		if(Noffs>0 & length(males>0)){
			#Determine the father
			fat<-sample(males,1)
			for(j in 1:Noffs){
				#cat(j,"\n")
				# Create the offspring
				pop<-c(pop,new("Leprechaun",parent1=i,parent2=fat))
			}
		}		
	}
  
	if(from!=CID){
	ALIVE<-c(ALIVE,(from):(CID-1))
	}
	### Everything should be written to a dataframe, to make sure we have all the values for ever and ever
  for(i in ALIVE){
    cat("\n",YR,"\t",pop[[i]]@ID,"\t",pop[[i]]@size,"\t",pop[[i]]@bvs,"\t",pop[[i]]@camemberts,"\t",pop[[i]]@sex,"\t",pop[[i]]@ARS,"\t",pop[[i]]@age,"\t",pop[[i]]@pID[1],"\t",pop[[i]]@pID[2],"\t",1,file=filename,append=TRUE)
  }
  for(i in DEAD){
    cat("\n",YR,"\t",pop[[i]]@ID,"\t",pop[[i]]@size,"\t",pop[[i]]@bvs,"\t",pop[[i]]@camemberts,"\t",pop[[i]]@sex,"\t",pop[[i]]@ARS,"\t",pop[[i]]@age,"\t",pop[[i]]@pID[1],"\t",pop[[i]]@pID[2],"\t",0,file=filename,append=TRUE)
  }
}

