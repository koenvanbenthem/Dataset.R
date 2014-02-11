#####################################################
### Overview of what we want to do/include/model ####
############## in our simulated dataset. ############

#####################################################
########## (1) 'Settings': ##########################
#####################################################
### Every individual is characterized by:
# Sex (constant through life) (M/F)
# Age/stage (changes through life) (x)
# Heritable phenotypic trait value (e.g. birth weight) (z) (changes/constant through life depending on trait)
# (Possibly: maternal effect m)
# (Possibly: breeding value a)

### Define relations between x, z and demographic processes (Possibly population density d included in functions) (for now, no sex differences)
# Survival(x,z,d) (logistic function)
# Growth (x,z,d) (either transition probability to next stage or absolute growth)
# Reproduction probability: p_repr(x,z,d) (logistic function)
# Number of offspring: n_offspring(x,z,d) (poisson distribution)
# Offspring size x distribution: x_offspring(x_mother,z_mother,x_father,z_father) (not required if working with stages) (prability density function)
# Offspring trait z distribution: z_offspring(x_mother,z_mother,x_father,z_father,a_mother,a_father,m) + rnorm(0,V_E) (prability density function)

### Environmental aspects
# Extra variation in z due to unexplained factors (i.e. V_E)) (assumed constant over time and constant over all individuals)
# Changes in the environment affecting survival(x,z), growth(x,z), p_repr(x,z), n_offspring(x,z), f_offspring(x_mother,z_mother)
#(i.e. changing selection)

#####################################################
########## (2) Define ancenstral population: ########
#####################################################
# n start individuals
# Start trait z distribution
# Start age/stage distribution
# Assign sexes to individuals
# (Possibly: start a and m values)

#####################################################
########## (3) In each time step: ###################
#####################################################
######## In following order:
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
