## Setting the Genetic maps

SDZ<-1
SDH<-1

############### Genetic determinism Z ################
dominance<-1 # for additive effects only, must be 0
overdominance<-0 # non-null values generate overdominance

nbLoci<-10 #number of loci controling the trait phenotype
nbAlleles<-10 #number of existing alleles per loci

gvaluesZ<-array(data=NA,dim=c(nbAlleles,nbAlleles,nbLoci),dimnames=list(paste("A",1: nbAlleles,sep=""),paste("A",1: nbAlleles,sep=""),paste("L",1:nbLoci,sep=""))) # Initialising a matrix that will contain the genotypic effects on the/a trait

for(L in 1:nbLoci)
{
  # Setting the effects for the homozygotes [all loci]
  effect<-abs(rnorm(n=1,mean=0,sd=SDZ))# alter the locus importance in a realistic way (many small-effect loci, few major loci)
  diag(gvaluesZ[,,L])<-2*rnorm(n=dim(gvaluesZ)[1],mean=0,sd=effect)
  # Setting the effects for the heterozygotes
  for(A in 1:(nbAlleles-1))# loop for off-diagonal = heterozygotes (additive and dominance effects)
  {
    for (D in (A+1):nbAlleles)
    {
      d<-dominance*runif(n=1,min=-0.5-overdominance,max=0.5+overdominance)
      gvaluesZ[A,D,L]<-(0.5-d)*gvaluesZ[A,A,L]+(0.5+d)*gvaluesZ[D,D,L] # mean of additive effects + dominance, over diagonal
      gvaluesZ[D,A,L]<-(0.5-d)*gvaluesZ[A,A,L]+(0.5+d)*gvaluesZ[D,D,L] # the same below diagonal    
    }
  }
}
############### Genetic determinism H ################
dominance<-1 # for additive effects only, must be 0
overdominance<-0 # non-null values generate overdominance

nbLoci<-10 #number of loci controling the trait phenotype
nbAlleles<-10 #number of existing alleles per loci

gvaluesH<-array(data=NA,dim=c(nbAlleles,nbAlleles,nbLoci),dimnames=list(paste("A",1: nbAlleles,sep=""),paste("A",1: nbAlleles,sep=""),paste("L",1:nbLoci,sep=""))) # Initialising a matrix that will contain the genotypic effects on the/a trait

for(L in 1:nbLoci)
{
  # Setting the effects for the homozygotes [all loci]
  effect<-abs(rnorm(n=1,mean=0,sd=SDH))# alter the locus importance in a realistic way (many small-effect loci, few major loci)
  diag(gvaluesH[,,L])<-2*rnorm(n=dim(gvaluesH)[1],mean=0,sd=effect)
  # Setting the effects for the heterozygotes
  for(A in 1:(nbAlleles-1))# loop for off-diagonal = heterozygotes (additive and dominance effects)
  {
    for (D in (A+1):nbAlleles)
    {
      d<-dominance*runif(n=1,min=-0.5-overdominance,max=0.5+overdominance)
      gvaluesH[A,D,L]<-(0.5-d)*gvaluesH[A,A,L]+(0.5+d)*gvaluesH[D,D,L] # mean of additive effects + dominance, over diagonal
      gvaluesH[D,A,L]<-(0.5-d)*gvaluesH[A,A,L]+(0.5+d)*gvaluesH[D,D,L] # the same below diagonal    
    }
  }
}
