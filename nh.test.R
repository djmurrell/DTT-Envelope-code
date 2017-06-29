install.packages("devtools")
library(devtools)
install_github('myllym/spptest', ref = 'no_fastdepth')

library(devtools)
library(spptest)
library(geiger)
library(phytools)
#library(doMC)
#help(spptest)

#registerDoMC()

###Number of species
spp_num<-160
realisations<-100
nsims<-1000

setwd("~/Dropbox/DTT")

spp_num<-60

burst<- -2

for(burst in c(-2, -4, -8))
for(spp_num in c(10,20,30,40,50,60,70,80,100,120,140,160,180,200))
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
				
	tree1<-pbtree(n=spp_num,scale=1,b=0.1, d=0.0, nsim=1)	
	x1<-fastBM(rescale(tree1,"EB",a=burst))
	#x1<-fastBM(tree1)
	d1<-dtt(tree1, x1, nsim=nsims, plot=T, calculateMDIp=T)

	if(d1$MDIp<0.05)
		mdi<-mdi+1
		
	cat("MDI-pvalue= ", d1$MDIp, "\n")
	cat("MDI = ", mdi/i, "\t")
	
	nh1<-nh.test(tree1, x1, regression.type="lm")
	
	if(nh1$coefficients[2,4]<0.05)
		nh<-nh+1
		
	cat("NH = ", nh/i, "\t")
	
###From Swensson book
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

	#But a specific test at one point has fp rate of ~5%
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
	
	#if(spp_num>11){
		#fp<-read.table(file="rank.envelope.tp.txt")
		#plot(fp$V1, fp$V2, ylim=c(0, 1), xlab="Species number", ylab="False positive rate")
		#lines(fp$V1, fp$V3, lwd=2)
		#y1<-fp$V1*0 +0.05
		#lines(fp$V1, y1, lty=2)
		#}

	
	if(res$p<0.05)
		num<-num+1
	

	#cat("MDI p-value = ", d1$MDIpVal, "\t")
	#cat("Rank p-value ", res$p, "\t")	
	cat("Rank test: ", num/i, "\n")	
		cat(100*i/realisations, "% complete\n\n")
	}
		
df<-data.frame(spp<-spp_num, mdi<-mdi/i, nh<-nh/i, pwise<-hyp/i, rank<-num/i, sims<-nsims, mc<-realisations, burst<-burst)

write.table(df, file="rank.envelope.nh.txt", append=TRUE, col.names=F, sep="\t", row.names=F)

}