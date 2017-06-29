#This is a modification of the MDI function in Geiger
#Here the MDI value for the data is compared to the MDI value for each of the simulations of the null model
#The data is 'significant' if it falls outside of either of the confidence intervals (ie is less than the 2.5% quantile
#or greater than the 97.5% quantile).
#
#
#This function is called by dtt1 (which is turn a modification from Geiger).
#
#


getMDIp1<-function (dttRes) 
{
    foo <- function(x) {
    	avedtt<-apply(dttRes$sim,1,mean)
        return(geiger:::.area.between.curves(x = dttRes$times, f1 = x, 
            f2 = avedtt))
    }
    avedtt<-apply(dttRes$sim,1,mean)
    mdis <- apply(dttRes$sim, 2, foo)
    
    mdis<-sort(mdis)
    #cat(mdis)
    
    pVal<-1
        	
    	cat(quantile(mdis, probs=c(0.025, 0.975)), '\n')
    	cat(dttRes$MDI, '\n')
    	
    #if(dttRes$MDI<quantile(mdis, probs=c(0.025, 0.975))[1])
    	#pVal<-0
    #if(dttRes$MDI>quantile(mdis, probs=c(0.025, 0.975))[2])
    	#pVal<-0
    	
    	
    	if(dttRes$MDI<mdis[0.025*(nsims)])
    		pVal<-0
    	else if(dttRes$MDI>mdis[0.975*(nsims)])
    		pVal<-0
    		
        return(pVal)
}
