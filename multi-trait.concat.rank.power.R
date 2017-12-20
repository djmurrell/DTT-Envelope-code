#
#
# Here we test the power to detect non-random trait evolution when we have two traits using the concatenated rank envelope test method.
# Two traits are given independent evolution rates, and we compute separate DTT curves. These curves are then concatenated and the rank test is performed on the single (concatenated) curve
# Test of pwer to correctly reject the null model is investigated by allowing one trait to have different rates of evolution to the (fixed rate) first trait

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
hyp<-0
mdi<-0
nh<-0

par(mfrow=c(3,1))

start<-1

	
for(i in start:realisations)
	{
				
	#Trees are simple Yule process models
	tree1<-pbtree(n=spp_num,scale=1,b=0.1, d=0.0, nsim=1)	
	
	#But we can rescale trait evolution to model early/late bursts
	x1<-fastBM(rescale(tree1,"EB",a=-4*log(2)))
	x2<-fastBM(rescale(tree1,"EB",a=burst))
	
	X1<-data.frame(x1=x1,x2=x2)
	#x1<-fastBM(tree1)
	
	#Call dtt for the classical DTT analyses
	d1.1<-dtt1(tree1, x1, nsim=nsims, plot=T, calculateMDIp=T)
	
	d1.2<-dtt1(tree1, x2, nsim=nsims, plot=T, calculateMDIp=T)

	dtt.sims<-rbind(d1.1$sim, d1.2$sim)
	dtt.obs<-rbind(d1.1$dtt, d1.2$dtt)
	dtt.obs<-c(d1.1$dtt, d1.2$dtt)
			
	sims<-dtt.sims
	sims<-as.matrix(sims)

	s1<-sims[-c(1),]
	
	r<-c(as.vector(d1.1$times), as.vector(d1.2$times+1))

	r<-as.vector(r[-c(1)])

	obs<-as.vector(dtt.obs)
	obs<-obs[-c(1)]

	c1<-list(r,obs, s1)
	names(c1)=c("r","obs","sim_m") 
	c2<-create_curve_set(c1)

	res<-rank_envelope(c2, alternative="less")
	ylim<-par("yaxp")
	r1<-res
	#####Plot rank envelope####
	plot(c(0,2), ylim[1:2], yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main="")
	x<-r1$r
	y1<-r1$upper
	y2<-r1$lower
	polygon(c(x, rev(x)), c(y1, rev(y2)), col="grey60", border=NA)
	lines(x, r1$data_curve, lwd=1)
	lines(x, r1$central_curve, lty=2)
	
	rl_y<-subset(r1$data_curve, r1$upper<r1$data_curve)
	
	rl_x<-subset(r1$r, r1$upper<r1$data_curve)
	
	points(rl_x, rl_y, cex=0.5)
	
	rl_y<-subset(r1$data_curve, r1$lower>r1$data_curve)
	
	rl_x<-subset(r1$r, r1$lower>r1$data_curve)
	
	points(rl_x, rl_y)

		
	if(res$p<0.05)
		num<-num+1
		
	cat("Rank test: ", num/i, "\n")	
		cat(100*i/realisations, "% complete\n\n")
		
		df<-data.frame(spp<-spp_num, rank<-num, sims<-nsims, mc<-i, burst<-burst)
		write.table(df, file="revision.tp.multivar.rank1.txt", append=TRUE, col.names=F, sep="\t", row.names=F)
	}
		
df<-data.frame(spp<-spp_num, rank<-num/i, sims<-nsims, mc<-realisations, burst<-burst)

write.table(df, file="revision.tp.multitvar.rank.300.txt", append=TRUE, col.names=F, sep="\t", row.names=F)

}
