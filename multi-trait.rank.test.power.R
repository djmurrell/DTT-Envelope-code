#
#	This code tests false and true positive rates for the MDI and rank envelope tests when there are two traits
#	In this case the relative disparity is computed using the modified dtt function of geiger, meaning the disparity value
#	is multi-varaite (averaging over both traits at once)
#
#	The rate of evolution for one trait is held constant with a half life of -4, and the half life for the second trait is allowed to vary
#
install.packages("devtools")
library(devtools)
install_github('myllym/spptest', ref = 'no_fastdepth')

library(devtools)
library(spptest)
library(geiger)
library(phytools)


#Set the number of simulated trees to test
realisations<-300

#Set the number of MC trait evolution simulations to run on each tree
nsims<-2500

################
#NEED to use CORRECTED source code from https://github.com/mwpennell/geiger-v2/blob/master/R/disparity.R
################
source("https://raw.githubusercontent.com/mwpennell/geiger-v2/master/R/disparity.R")

#Modified dtt function to call modified MDI function and allow user to change y-axis
source("https://raw.githubusercontent.com/djmurrell/DTT-Envelope-code/master/dtt1.R")

#Modified code for two sided MDI test
source("https://raw.githubusercontent.com/djmurrell/DTT-Envelope-code/master/getMDI1.R")

source("https://raw.githubusercontent.com/djmurrell/DTT-Envelope-code/master/rank_dtt.R")


#variable for number of species at the tips
spp_num<-10

#Determine the rate of change of trait evolution: 
#		half_life<0 for a decelerating rate
#		half_life>0 for an accelerating rate
#		half_life=0 for Brownian evolution (ie for estimating false positive rates)
#
#	Changing this variable only changesthe half life for the second trait 

half_life<- 0

#We can loop through all species numbers
#for(spp_num in c(10,20,30,40,50,60,70,80,100))

#And/or we could loop through the half lives and doubling times to be investigated
#for(half_life in c(-4, -2, 0))

{
###Number of datasets to simulate and test against null model

burst<-half_life*log(2)
cat(burst)
num<-0
pwise<-0
mdi<-0


par(mfrow=c(3,1))

start<-1

	
for(i in start:realisations)
	{
				
	#Trees are simple Yule process models
	tree1<-pbtree(n=spp_num,scale=1,b=0.1, d=0.0, nsim=1)	
	
	#But we can rescale trait evolution to model early/late bursts
	#SIMULATE EVOLUTION OF TWO TRAITS
	x1<-fastBM(rescale(tree1,"EB",a=-4*log(2)))
	x2<-fastBM(rescale(tree1,"EB",a=burst))
	
	X1<-data.frame(x1=x1,x2=x2)
	#x1<-fastBM(tree1)
	
	#Call dtt for the classical DTT analyses
	d1<-dtt1(tree1, X1, nsim=nsims, plot=T, calculateMDIp=T)

##Note we need to decide if we are using the 1-tailed (p<0.05) or 2-tailed (p<0.025) version
	if(d1$MDIp<0.05)
		mdi<-mdi+1
		
	cat("MDI-pvalue= ", d1$MDIp, "\n")
	cat("MDI = ", mdi/i, "\t")
	

#####################################################################
#
#		Call function to perform rank envelope test
#
#####################################################################

	res<-rank_env_dtt(d1, test="less")

		
	if(res$p<0.05)
		num<-num+1
		
	cat("Rank test: ", num/i, "\n")	
		cat(100*i/realisations, "% complete\n\n")
		
		df<-data.frame(spp<-spp_num, mdi<-mdi, rank<-num, sims<-nsims, mc<-i, burst<-burst)
		write.table(df, file="revision.tp.multitrait.txt", append=TRUE, col.names=F, sep="\t", row.names=F)
	}
		
df<-data.frame(spp<-spp_num, mdi<-mdi/i, rank<-num/i, sims<-nsims, mc<-realisations, burst<-burst)

write.table(df, file="revision.tp.multitrait.300.txt", append=TRUE, col.names=F, sep="\t", row.names=F)

}
