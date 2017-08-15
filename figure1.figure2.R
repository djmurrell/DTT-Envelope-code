install.packages("devtools")
library(devtools)
install_github('myllym/spptest', ref = 'no_fastdepth')

library(devtools)
library(spptest)
library(geiger)
library(phytools)


#Set the number of simulated trees to test
realisations<-100

#Set the number of MC trait evolution simulations to run on each tree
nsims<-5000

setwd("~/Dropbox/DTT")

#variable for number of species at the tips
spp_num<-20

#When burst > 0 we model a late burst in morphological evolution; 
#When burst < 0 we model an early burst in trait evolution
#When burst = 0 we model Brownian evolution
burst<- 0

#We can loop through all species numbers

#for(spp_num in c(10,20,30,40,50,60,70,80,100,120,140,160,180,200))
{
###Number of datasets to simulate and test against null model

num<-0
pwise<-0
hyp<-0
mdi<-0
nh<-0

par(mfrow=c(3,1))

for(i in 1:realisations)
	{
				
	#Trees are simple Yule process models
	tree1<-pbtree(n=spp_num,scale=1,b=0.1, d=0.0, nsim=1)	
	
	#But we can rescale trait evolution to model early/late bursts
	x1<-fastBM(rescale(tree1,"EB",a=burst))
	#x1<-fastBM(tree1)
	
	#Call dtt for the classical DTT analyses
	d1<-dtt(tree1, x1, nsim=nsims, plot=T, calculateMDIp=T)

	if(d1$MDIp<0.05)
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

	#Point wise (multiple) testing has fp rate of ~30%
	cat("pwise = ", pwise/i, "\t")

	#But a specific test at one point/node has fp rate of ~5%
	if( sum(ifelse(t1[seq(2,spp_num-1,1)]>0.975*(nsims+1), 1, 0)) >0)
		hyp=hyp+1
	else if( sum(ifelse(t1[seq(2,spp_num-1,1)]<0.025*(nsims+1), 1, 0)) >0)
		hyp=hyp+1
	cat("hypothesis = ", hyp/i, "\t")

#####################################################################

	sims<-d1$sim
	sims<-as.matrix(sims)

	s1<-sims[-c(1, spp_num),]
	
	r<-d1$times[-c(1,(spp_num))]
	r<-as.vector(r)

	obs<-as.vector(d1$dtt)
	obs<-obs[-c(1,(spp_num))]

	c1<-list(r,obs, s1)
	names(c1)=c("r","obs","sim_m") 
	c2<-create_curve_set(c1)

	res<-rank_envelope(c2)

	#####Plot rank envelope####
	plot(res)
		
	if(res$p<0.05)
		num<-num+1
		
	cat("Rank test: ", num/i, "\n")	
		cat(100*i/realisations, "% complete\n\n")
	}
		
df<-data.frame(spp<-spp_num, mdi<-mdi/i, nh<-nh/i, pwise<-hyp/i, rank<-num/i, sims<-nsims, mc<-realisations, burst<-burst)

write.table(df, file="rank.envelope.nh.txt", append=TRUE, col.names=F, sep="\t", row.names=F)

}