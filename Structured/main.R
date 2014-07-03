##### In this file we will try to come up with some things ####
###############################################################

########################
####### Type of program
## Written in R
## Object Oriented
## Uses parallel computations with library("Rmpi")
#######################

#######################
# Main object "Leprechaun": an individual of our species of interest, the Irish Leprechaun, small in size, but with great powers. To make it more French for Timothee, we assume that their favourite food is camembert. Leprechaun is not very choosy, and mates completely random.
##########################################

# The object contains the following values

# Static [these numbers do not change after initialisation]
#  V ID (unique identifier number) [integer]
#  V pID (two numbers referring to the parents of the Leprechaun, if none: NA) [vector of two integers]
#  V Year of birth (timestep at which the individual was born)	[integer]
#  V Genome (?) (two vectors of length N coding for both chromosomes of N loci in the genome.) [two vectors of N integers]
#  V Heritable phenotypic trait value of interest (z) (e.g. birth weight) (changes/constant through life depending on trait)
#  - (Possibly: breeding value A))
#  - Rather simulate physically independent loci, otherwise we need to simulate recombination on the chromosomes. Over short time periods, 100 recombinations fragments (ie independent loci) sounds realistic.
#  V Explicit coding of traits by many independent diploid loci. The simplest model: z = mean + sum_loci(a1_locus + a2_locus) + environment. One can add explicit dominance and epistasis, as well as interactions with environment. 
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

############### Setting folder for storage ###############
wd<-"/Users/koen/Dataset.R/Structured"
setwd(wd)
folder<-"Test3"

converter<-paste("Data/",folder,"/conv.csv",sep="")
dir.create(file.path("Data", folder))


################ Setting parameters ########################
norep<-1 # Number of replicates
SS1<-seq(1,5,1) # Linear survival selection
SS2<-c(-0.01) # Quadratic survival selection
fs1<-c(0.1) # Linear fertility selection
fs2<-c(-0.01) # Quadratic fertility selection
cs<-c(0.1) # Camembert selection
sdc<-c(1000) #SDcamembert

############### Creating all tasks #########################
tasks<-expand.grid(rep=1:norep,SS1=SS1,SS2=SS2,fs1=fs1,fs2=fs2,cs=cs,sdc=sdc)
tasks$seed<-sample(1:length(tasks$rep))
tasks$genseed<-rep(789,length(tasks$rep))
tasks$filename<-paste("Data/",folder,"/pop_seed_",tasks$seed,"genseed_",tasks$genseed,"SS1_",tasks$SS1*10,"_SS2_",tasks$SS2*100,"_FS1_",tasks$fs1*10,"_FS2_",tasks$fs2*100,"_CAMSEL_",tasks$cs*10,"_SDCAM_",tasks$sdc,".csv",sep="")

write.csv(tasks,file=converter)

############## Starting the slaves ######################
library("Rmpi")


############## Initialising the slaves by learning them the classes ####################

############## Performing the calculations ###################################
# Notice we just say "give us all the slaves you've got."
mpi.spawn.Rslaves()

if (mpi.comm.size() < 2) {
  print("More slave processes are required.")
  mpi.quit()
}

.Last <- function(){
  if (is.loaded("mpi_initialize")){
    if (mpi.comm.size(1) > 0){
      print("Please use mpi.close.Rslaves() to close slaves.")
      mpi.close.Rslaves()
    }
    print("Please use mpi.quit() to quit R")
    .Call("mpi_finalize")
  }
}

# Function the slaves will call to perform
slave <- function() {
  source("Leprechaun.R",local=TRUE)
  # Note the use of the tag for sent messages: 
  #     1=ready_for_task, 2=done_task, 3=exiting 
  # Note the use of the tag for received messages: 
  #     1=task, 2=done_tasks 
  junk <- 0 
  
  done <- 0 
  while (done != 1) {
    # Signal being ready to receive a new task 
    mpi.send.Robj(junk,0,1) 
    
    # Receive a task 
    task <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag()) 
    task_info <- mpi.get.sourcetag() 
    tag <- task_info[2]
    
    if (tag == 1) {
      set.seed(task$genseed)
      source("Genes.R",local=TRUE)
      set.seed(task$seed)
      source("Time.R",local=TRUE)
      
    } else if (tag == 2) {
      done <- 1
    }
    # We'll just ignore any unknown messages
  }
  
  mpi.send.Robj(junk,0,3)
}

# Send the function to the slaves
mpi.bcast.Robj2slave(slave)

# Call the function in all the slaves to get them ready to
# undertake tasks
mpi.bcast.cmd(setwd(wd))
mpi.bcast.cmd(slave())

# Teach them what a leprechaun is
#mpi.bcast.cmd(source("Leprechaun.R",local=TRUE))

junk <- 0 
closed_slaves <- 0 
n_slaves <- mpi.comm.size()-1 

#convert the tasks data frame into a list
tasks<-split(tasks, rownames(tasks))

while (closed_slaves < n_slaves) { 
  # Receive a message from a slave 
  message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag()) 
  message_info <- mpi.get.sourcetag() 
  slave_id <- message_info[1] 
  tag <- message_info[2] 
  
  if (tag == 1) { 
    # slave is ready for a task. Give it the next task, or tell it tasks 
    # are done if there are none. 
    if (length(tasks) > 0) { 
      # Send a task, and then remove it from the task list 
      mpi.send.Robj(tasks[[1]], slave_id, 1); 
      tasks[[1]] <- NULL 
    } 
    else { 
      mpi.send.Robj(junk, slave_id, 2) 
    } 
  } 
  else if (tag == 2) { 
    # The message contains results. Do something with the results. 
    # Store them in the data structure -- no, was already done before

  } 
  else if (tag == 3) { 
    # A slave has closed down. 
    closed_slaves <- closed_slaves + 1 
  } 
} 



mpi.close.Rslaves()
mpi.quit(save="no")