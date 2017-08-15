install.packages("devtools")
library(devtools)
install_github('myllym/spptest', ref = 'no_fastdepth')

library(spptest)
library(geiger)

#Number of simulations to sample from null model of Brownian evolution.
#Recommendation for at least 2500

nsims=5000

#Slight alteration of dtt function in geiger
source("dtt1.R")


####################################################
#
#	Generic function to compute rank envelope test
#
####################################################

rank_env_dtt<-function(x, Plot=T)
{
		
	spp_num<-length(x$times)		
	sims<-x$sim
	sims<-as.matrix(sims)

	s1<-sims[-c(1),]
	
	r<-x$times[-c(1)]

	r<-as.vector(r)

	obs<-as.vector(x$dtt)
	obs<-obs[-c(1)]

	c1<-list(r,obs, s1)
	names(c1)=c("r","obs","sim_m") 
	c2<-create_curve_set(c1)

	res<-rank_envelope(c2)
	
	if(Plot==T)
	plot(res, xlab="Relative time", ylab="Disparity", main="")	
	return(res)	
}

####################################################
#
#	Re-test of cichlid morphological data from Feilich (2017)
#   See http://onlinelibrary.wiley.com/doi/10.1111/evo.13021/abstract
#
####################################################

#Data file downloaded from http://datadryad.org/resource/doi:10.5061/dryad.h4k6f/7
load("06.23.2016_DTT_wkspc.RData")


pdf("comparison_data.supp.pdf")

	layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE), widths=c(2,2,2),heights=c(2,2,2))
	
	#Plot pointwise DTT then rank envelope DTT for body shape of females
	rn<-row.names(bodydata_F)
	t1<-data.frame(t1=bodydata_F, row.names=rn)
	d1<-dtt1(bodytree_F, bodydata_F, nsim=nsims, Ylim=c(0,4))
	
	ylim<-par("yaxp")
	
	#Compute the rank envelope
	r1<-rank_env_dtt(d1, Plot=F)

	#Note that the p-value is a range because the ranks will almost always lead to some ties
	r1$p_interval

	plot(c(0,1), ylim[1:2], yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main="")
	x<-r1$r
	y1<-r1$upper
	y2<-r1$lower
	polygon(c(x, rev(x)), c(y1, rev(y2)), col="grey60", border=NA)
	lines(x, r1$data_curve, lwd=1)
	lines(x, r1$central_curve, lty=2)
	
	#########################################

	#Plot pointwise DTT then rank envelope DTT for caudal fin shape of females	
	rn<-row.names(caudaldata_F)
	t1<-data.frame(t1=caudaldata_F, row.names=rn)
	d1<-dtt1(caudaltree_F, caudaldata_F, nsim=nsims, Ylim=c(0,4))

	ylim<-par("yaxp")
	
	#Compute the rank envelope
	r1<-rank_env_dtt(d1, Plot=F)

	#Note that the p-value is a range because the ranks will almost always lead to some ties
	r1$p_interval

	plot(c(0,1), ylim[1:2], yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main="")
	x<-r1$r
	y1<-r1$upper
	y2<-r1$lower
	polygon(c(x, rev(x)), c(y1, rev(y2)), col="grey60", border=NA)
	lines(x, r1$data_curve, lwd=1)
	lines(x, r1$central_curve, lty=2)
	
	#########################################
	
	#Plot pointwise DTT then rank envelope DTT for dorsal fin shape of females
	rn<-row.names(dorsaldata_F)
	t1<-data.frame(t1=dorsaldata_F, row.names=rn)
	d1<-dtt1(dorsaltree_F, dorsaldata_F, nsim=nsims, Ylim=c(0,4))
	
	ylim<-par("yaxp")
	
	#Compute the rank envelope
	r1<-rank_env_dtt(d1, Plot=F)

	#Note that the p-value is a range because the ranks will almost always lead to some ties
	r1$p_interval

	plot(c(0,1), ylim[1:2], yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main="")
	x<-r1$r
	y1<-r1$upper
	y2<-r1$lower
	polygon(c(x, rev(x)), c(y1, rev(y2)), col="grey60", border=NA)
	lines(x, r1$data_curve, lwd=1)
	lines(x, r1$central_curve, lty=2)
	
dev.off()