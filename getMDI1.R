#This is a modification of the MDI function in Geiger
#Here the p-value is computed as in Slater GJ, Price, SA, Santini, F, Alfaro, MA. 2010. Proceedings of the Royal Society B. 277: 3097 -3104.
#However, the original code is for a one-tailed test, and the modification below is for a two-tailed test.
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
