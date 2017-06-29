library(devtools)
install_github('myllym/spptest', ref = 'no_fastdepth')

library(devtools)
library(spptest)
library(geiger)
library(phytools)

spp_num<-100

nsims<-5000

setwd("~/Dropbox/DTT")



pdf("example.dtt.pdf")
	#par(mfrow=c(3,2))
	layout(matrix(c(1,6,2,5,3,4), 3, 2, byrow = TRUE))

source("dtt1.R")
for(burst in c(-8,-4,-2,2,4,8))
	{
		tree1<-pbtree(n=spp_num,scale=1,d=0.0, nsim=1)	
		x1<-fastBM(rescale(tree1,"EB",a=burst))
		d1<-dtt1(tree1, x1, nsim=nsims, plot=T)
		
		mylabel = bquote(a == .(format(burst, digits = 1)))

		text(0.8, 1.5, labels=mylabel)
		
	}
	
dev.off()	