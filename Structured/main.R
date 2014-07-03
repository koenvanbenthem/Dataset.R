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

# Loading the class and the genetic map
source("Leprechaun.R")
source("Genes.R")

############### Selection parameters ###############
folder<-"StrongerSSMoreAndAlsoASecondOrderTerm2"
converter<-paste("Data/",folder,"/conv.csv",sep="")
dir.create(file.path("Data", folder))
cat("filename\tSS1\tSS2\tFS1\tFS2\tCAMSEL\tSDCAM",file=converter,append=FALSE)

for(SurvivalSelection1 in seq(1,5,1)){
for(SurvivalSelection2 in c(-0.01)){
for(fertilitySelection1 in c(0.1)){
for(fertilitySelection2 in c(-0.01)){
for(camembertSelection in c(0.1)){
for(SDCamembert in c(1000)){
  
#SurvivalSelection1<-0.1 #linear coefficient on a logit scale for Survival ~ ... + size +size^2
#SurvivalSelection2<-(-0.01) #quadratic coefficient on a logit scale for Survival ~ ... + size + size^2; negative value=balancing selection

#fertilitySelection1<-0.1 #linear coefficient on a log scale for reproduction ~ ... + size + size^2
#fertilitySelection2<-(-0.01) #quadratic coefficient on a log scale for reproduction ~ ... + size + size^2; negative value=balancing selection

#camembertSelection<-0.1
survivalPenaltyForRepro<-0

SDZ<-1
SDH<-1
####################################################

################ Environmental parameters ##########
MeanCamembert<-5000
#SDCamembert<-1000

lowBoundGrowth<-0.99 # minimal growth rate
highBoundGrowth<-1.1 # maximal growth rate

filename<-paste("Data/",folder,"/pop_SS1_",SurvivalSelection1*10,"_SS2_",SurvivalSelection2*100,"_FS1_",fertilitySelection1*10,
                "_FS2_",fertilitySelection2*100,"_CAMSEL_",camembertSelection*10,"_SDCAM_",SDCamembert,".csv",sep="")
cat("\n",filename,"\t",SurvivalSelection1,"\t",SurvivalSelection2,"\t",fertilitySelection1,"\t",fertilitySelection2,"\t",camembertSelection,"\t",SDCamembert,file=converter,append=TRUE)



################ Random seed #######################
set.seed(12)

################ Counter for the IDs ###############
CID<-as.integer(1) 


################ Counter for the current year ######
YR<-0 

################ Base population parameters ########

MeanBirthSize<-10
# needed to make survival ~ size relative on age
meanGrowth<-prod(runif(10000,lowBoundGrowth,highBoundGrowth))^(1/10000)
MeanRepro<-2

source("Time.R")

}}}}}}
