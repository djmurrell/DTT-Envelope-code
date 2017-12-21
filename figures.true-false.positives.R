install.packages("devtools")
library(devtools)
install_github('myllym/spptest', ref = 'no_fastdepth')

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
#		half_life=0 for Brownian evolution

half_life<- 0

#We can loop through all species numbers
#for(spp_num in c(10,20,30,40,50,60,70,80,100,120,140,160,180,200))

#And/or we could loop through the half lives and doubling times to be investigated
#for(half_life in c(-4, -5, -6, -5, 4, 5, 6))

{
###Number of datasets to simulate and test against null model

burst<-half_life*log(2)
cat(burst)
num<-0
pwise<-0
mdi<-0
nh<-0

par(mfrow=c(3,1))

start<-1

	
for(i in start:realisations)
	{
				
	#Trees are simple Yule process models
	tree1<-pbtree(n=spp_num,scale=1,b=0.1, d=0.0, nsim=1)	
	
	#But we can rescale trait evolution to model early/late bursts
	x1<-fastBM(rescale(tree1,"EB",a=burst))
	#x1<-fastBM(tree1)
	
	#Call dtt for the classical DTT analyses
	d1<-dtt1(tree1, x1, nsim=nsims, plot=T, calculateMDIp=T)

##Note we need to decide if we are using the 1-tailed (p<0.05) or 2-tailed (p<0.025) version
	if(d1$MDIp<0.025)
		mdi<-mdi+1
		
	cat("MDI-pvalue= ", d1$MDIp, "\n")
	cat("MDI = ", mdi/i, "\t")
	
	nh1<-nh.test(tree1, x1, regression.type="lm")
	
	if(nh1$coefficients[2,4]<0.05)
		nh<-nh+1
		
	cat("NH = ", nh/i, "\t")
	
###This is the pointwise testing approach but suffers
###from problems with multiple testing as it is not a global test.
###However, it works very well for testing a specific hypothesis
###relating to a specific time period. The danger here is that
###the specific time period is determined post hoc, and a post hoc
###intepretation is given. 

 	t1<-apply(cbind(d1$dtt, d1$sim), MARGIN=1, rank)[1,]
	#t1<0.03*(nsims+1)
	#t1>0.97*(nsims+1)
	
	
	if( sum(ifelse(t1>0.975*(nsims+1), 1, 0)) >0)
		pwise=pwise+1
	else if( sum(ifelse(t1<0.025*(nsims+1), 1, 0)) >0)
		pwise=pwise+1

	#Point wise (multiple) testing 
	cat("pwise = ", pwise/i, "\t")

#####################################################################
#
#		Call function to perform rank envelope test
#
#####################################################################

	res<-rank_env_dtt(d1)

		
	if(res$p<0.05)
		num<-num+1
		
	cat("Rank test: ", num/i, "\n")	
		cat(100*i/realisations, "% complete\n\n")
		
		df<-data.frame(spp<-spp_num, mdi<-mdi, nh<-nh, pwise<-pwise, rank<-num, sims<-nsims, mc<-i, burst<-burst)
		
#Collect data after each test tree -this allows the user to check the true/false positive rates are settling down
		write.table(df, file="revision.figure1.2-tailed.p0.power.txt", append=TRUE, col.names=F, sep="\t", row.names=F)
	}
#Write the true/false positive rates after all realisations for each parameter set have been completed		
	df<-data.frame(spp<-spp_num, mdi<-mdi/i, nh<-nh/i, pwise<-pwise/i, rank<-num/i, sims<-nsims, mc<-realisations, burst<-burst)
	write.table(df, file="revision.figure1.2-tailed.300.txt", append=TRUE, col.names=F, sep="\t", row.names=F)

}
