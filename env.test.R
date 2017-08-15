
#################################################################################################

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
nsims<-1999

setwd("~/Dropbox/DTT")

for(spp_num in c(10,20,30,40,50,60,70,80,100,120,140,160,180,200))
{
###Number of datasets to simulate and test against null model


burst<- 4

num<-0
pwise<-0
hyp<-0
mdi<-0

par(mfrow=c(2,1))

for(i in 1:realisations)
	{
				
	tree1<-pbtree(n=spp_num,scale=1,d=0.0, nsim=1)	
	x1<-fastBM(rescale(tree1,"EB",a=burst))
	#x1<-fastBM(tree1)
	d1<-dtt(tree1, x1, nsim=nsims, plot=T, calculateMDIp=T)

	if(d1$MDIp<0.05)
		mdi<-mdi+1
		
	cat("MDI-pvalue= ", d1$MDIp, "\n")
	cat("MDI = ", mdi/i, "\t")
	
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
		
df<-data.frame(spp<-spp_num, mdi<-mdi/i, pwise<-hyp/i, rank<-num/i, sims<-nsims, mc<-realisations, burst<-burst)

write.table(df, file="rank.envelope.tp1.txt", append=TRUE, col.names=F, sep="\t", row.names=F)

}


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
	#s1<-sims
	
	r<-x$times[-c(1)]
	#r<-x$times
	r<-as.vector(r)

	obs<-as.vector(x$dtt)
	obs<-obs[-c(1)]
	#obs<-obs

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

nsims=5000
source("dtt1.R")
	data(geospiza)
	#par(mfrow=c(1,1))

	pdf("comparison_data3.pdf")
	#par(mar = c(0,0,0,0))
	#par(mfrow=c(2,2))
	layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE), widths=c(2,2,2),heights=c(2,2,2))
	#pdf("geospiza1.pdf")
	#d1<-dtt(geospiza$phy, geospiza$dat[,3], plot=T, nsim=nsims)	
	d1<-dtt1(geospiza$phy, geospiza$dat[,3], plot=T, nsim=nsims, Ylim=c(0,2))	
	#dev.off()
	
	ylim<-par("yaxp")
	r1<-rank_env_dtt(d1, Plot=F)

	#pdf("whales2.pdf")
	plot(c(0,1), ylim[1:2], yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main="")
	x<-r1$r
	y1<-r1$upper
	y2<-r1$lower
	polygon(c(x, rev(x)), c(y1, rev(y2)), col="grey60", border=NA)
	lines(x, r1$data_curve, lwd=1)
	lines(x, r1$central_curve, lty=2)
	#dev.off()

#########################################################################
#	
#Test whales data within geiger and from Slater et al. (2010)
#
#########################################################################
data(whales)
#par(mfrow=c(2,1))	
	#pdf("whales1.pdf")
	
	d1<-dtt1(whales$phy, whales$dat[,2], plot=T, nsim=nsims, calculateMDIp=T, Ylim=c(0,1))
	#dev.off()
	ylim<-par("yaxp")
	
	r1<-rank_env_dtt(d1, Plot=F)

	#pdf("whales2.pdf")
	plot(c(0,1), ylim[1:2], yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main="")
	x<-r1$r
	y1<-r1$upper
	y2<-r1$lower
	polygon(c(x, rev(x)), c(y1, rev(y2)), col="grey60", border=NA)
	lines(x, r1$data_curve, lwd=1)
	lines(x, r1$central_curve, lty=2)
	#dev.off()


#########################################################################
#		
#See https://github.com/simjoly/CourseComparativeMethods/blob/master/lecture7/Diversification.Rmd	
# Original paper: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0110618
#########################################################################

paramotree <-read.nexus("JamesoniaPartFind_MCC.tre")	
paramodata <-read.csv("JamesoniaTraits2014Analysis.csv")
# Remove species for which we don't have complete data
paramodata <- paramodata[!is.na(paramodata$Name_in_tree),]
# Remove species in the tree that are not in the data matrix
species.to.exclude <- paramotree$tip.label[!(paramotree$tip.label %in% 
                                                   paramodata$Name_in_tree)]
paramotree <- drop.tip(paramotree,species.to.exclude)
rm(species.to.exclude)


# Name the rows of the data.frame with the species codes used as tree labels
rownames(paramodata) <- paramodata$Name_in_tree
# Remove unecessary variables
paramodata <- paramodata[,-c(1:3,7:14)]
# Order the data in the same order as the tip.label of the tree. In the
# present example, this was already the case.
paramodata <- paramodata[paramotree$tip.label,]

# Replace NAs by the mean for the calculations
paramodata[is.na(paramodata[,"Altitude"]),"Altitude"] <- mean(paramodata[,"Altitude"],na.rm=TRUE)
# dtt plot
#nsims=10000
	#par(mfrow=c(1,1))
	
	#pdf("jamesonia1.pdf")
	#d1<- dtt(paramotree,paramodata['Altitude'], index="avg.sq", nsim=nsims,calculateMDIp=T)
    d1<- dtt1(paramotree,paramodata['Altitude'], index="avg.sq", nsim=nsims, calculateMDIp=T, Ylim=c(0,2))
    #dev.off()
    ylim<-par("yaxp")
    
    #pdf("test.pdf", paper="a4")
    #par(mfrow=c(2,2))
	r1<-rank_env_dtt(d1, Plot=F)
	#for(i in 1:4) plot(r1)
	#dev.off()
	
	#pdf("jamesonia2.pdf")
	plot(c(0,1), ylim[1:2], yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main="")
	x<-r1$r
	y1<-r1$upper
	y2<-r1$lower
	polygon(c(x, rev(x)), c(y1, rev(y2)), col="grey60", border=NA)
	lines(x, r1$data_curve, lwd=1)
	lines(x, r1$central_curve, lty=2)

dev.off()

  




#Test caudata datasetin geieger
data(caudata)
nsims=1000
par(mfrow=c(2,1))
	d1<-dtt1(caudata$phy, caudata$dat, plot=T, nsim=nsims)	
	r1<-rank_env_dtt(d1)		
	
		
#Test primate dataset in geieger
data(primates)
source("dtt1.R")
	par(mfrow=c(2,1))
	nsims=1000
	d1<-dtt1(primates$phy, primates$dat, plot=T, nsim=nsims)	
	r1<-rank_env_dtt(d1)
	
	
	data(caniformia)
	par(mfrow=c(2,1))
	nsims=1000
	d1<-dtt(caniformia$phy, caniformia$dat, plot=T, nsim=nsims)	
	r1<-rank_env_dtt(d1)	
	
	  
    
	
	m1<-read.tree("monkeys.rtf")
	m1_traits<-read.csv("monkey.traits.csv",row.names=1)
	x<-m1_traits
	#m1_traits<-m1_traits[,1-2]
	
	rn<-row.names(x)
	m2<-data.frame(pc=m1_traits[,1], row.names=rn)
	
	par(mfrow=c(2,1))
	d1<-dtt(m1, m2, plot=T, nsim=10000)
	r1<-rank_env_dtt(d1)
	

setwd("~/Dropbox/DTT")	
	g1<-read.tree("Geophagini_tree.phy")
	g1_data<-read.csv("LOGsl.csv", row.names=1)
	g1_data<-read.csv("LOGvariables.csv", row.names=1)
	g1_data<-as.matrix(g1_data)
	par(mfrow=c(2,1))
	nsims=1000	
		d1<-dtt(g1, g1_data[,1], plot=T, nsim=nsims)	
		r1<-rank_env_dtt(d1)
	
#Test turtles dataset in geieger
data(chelonia)
par(mfrow=c(2,1))
	nsims=1000
	d1<-dtt(chelonia$phy, chelonia$dat, plot=T, nsim=nsims)	
	r1<-rank_env_dtt(d1)
	
	
	
	
	
	spp_num<-length(chelonia$dat)		
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

	plot(res)

		