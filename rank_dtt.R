####################################################
#
#	Generic function to compute rank envelope test
#
#	two tailed test: test="two.sided"
#	one sided tests: test="less" OR test="greater"
#
####################################################

rank_env_dtt<-function(x, Plot=T, test="two.sided")
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

	res<-rank_envelope(c2, alternative=test)
	
	if(Plot==T)
	plot(res, xlab="Relative time", ylab="Disparity", main="")	
	return(res)	
}
