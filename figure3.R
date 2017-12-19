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




#########################################################################
#	
#				Geospiza in geiger		
#
#########################################################################


	data(geospiza)
	
	pdf("comparison_data.final.pdf")

	layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE), widths=c(2,2,2),heights=c(2,2,2))
	
	#Compute and plot the DTT with pointwise envelope
	d1<-dtt1(geospiza$phy, geospiza$dat[,3], plot=T, nsim=nsims, Ylim=c(0,2))	
	
	ylim<-par("yaxp")
	
	#Compute the rank envelope
	r1<-rank_env_dtt(d1, Plot=F)

	#Note that the p-value is a range because the ranks will almost always lead to some ties
	r1$p_interval


	#We'll use a bespoke plot to make it look 'pretty'
	plot(c(0,1), ylim[1:2], yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main="")
	x<-r1$r
	y1<-r1$upper
	y2<-r1$lower
	polygon(c(x, rev(x)), c(y1, rev(y2)), col="grey60", border=NA)
	lines(x, r1$data_curve, lwd=1)
	lines(x, r1$central_curve, lty=2)


#########################################################################
#	
#Test whales data within geiger and from Slater et al. (2010)
#
#########################################################################
data(whales)
	
	#Compute and plot the DTT with pointwise envelope
	d1<-dtt1(whales$phy, whales$dat[,2], plot=T, nsim=nsims, calculateMDIp=T, Ylim=c(0,1))

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


#########################################################################
#	
#Test cichlid data from Feilich (2016) on evolution of anal fin shape
#See http://onlinelibrary.wiley.com/doi/10.1111/evo.13021/abstract
#
#########################################################################

#Data file downloaded from http://datadryad.org/resource/doi:10.5061/dryad.h4k6f/7
load("06.23.2016_DTT_wkspc.RData")

	rn<-row.names(analdata_F)
	t1<-data.frame(t1=analdata_F[,1-7], row.names=rn)
	
	#Compute and plot the DTT with pointwise envelope
	d1<-dtt1(analtree_F, t1, plot=T, nsim=nsims, calculateMDIp=F, Ylim=c(0,6))

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
