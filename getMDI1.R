#This is a modification of the MDI function in Geiger
#Here the MDI value for the data is compared to the MDI value for each of the simulations of the null model
#The data is 'significant' if it falls outside of either of the confidence intervals (ie is less than the 2.5% quantile
#or greater than the 97.5% quantile).
#
#
#This function is called by dtt1 (which is turn a modification from Geiger).
#
#


getMDIp<-function(dttRes) {
	foo<-function(x) {
		return(.area.between.curves(x= dttRes$times, f1=x, f2=dttRes$dtt))
	}
	mdis<-apply(dttRes$sim,2,foo)
    
    #Two sided test
	p1<-length(which(mdis>=0))/length(mdis)
    p2<-length(which(mdis<=0))/length(mdis)
    
    pVal<-min(p1,p2)
	return(pVal)
}
