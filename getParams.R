## HMM (or CBS) distribution estimation function
## functionized 6/6/12
## added better alternative parameter estimation 7/18/12
## alyssa frazee

source("~/Documents/hopkins/research/_ebHMM-Rfunctions/locfdrFit.R") 
source("~/Documents/hopkins/research/_ebHMM-Rfunctions/paramHelpers.R")
getParams <- function(tstats, plots=FALSE, plotfile=NULL, verbose=F)
{
	if(plots & !is.null(plotfile)) pdf(file=plotfile)
	pi.z = sum(tstats==0)/length(tstats)
	nzt = tstats[tstats!=0]
	
	if(plots){
		hist(nzt,freq=FALSE,breaks=100,xlab="t",main="Histogram of nonzero t-statistics")
		zzaxis = seq(-5,5,length=1000)
		lines(zzaxis,dnorm(zzaxis),col="red")
		legend("topright","N(0,1)",lty=1,col="red",bty="n")
	}
	plotarg = ifelse(plots,1,0)
	fdrmodel.defaults = locfdrFit(nzt,plot=plotarg,verbose=F,main="initial locfdr")
	if(fdrmodel.defaults$needsfix==1) {
		fdrmodel <- locfdrFit(nzt,mlest=c(fdrmodel.defaults$mlest.lo, fdrmodel.defaults$mlest.hi),verbose=F,plot=plotarg,main="locfdr, incorporating mlest")
		if(fdrmodel$needsfix==1) stop("problem with mlest parameters in locfdr")
	}
	if(fdrmodel.defaults$needsfix==0) fdrmodel <- fdrmodel.defaults
	
	pi.0.tmp = fdrmodel$fp0[3,3]
	if(pi.0.tmp>=1) {
		warning("Proportion of null nonzero bases estimated to be greater than 1.  Defaulting to 0.999.")
		pi.0.tmp <- 0.999
	}
	pi.0 = pi.0.tmp*(1-pi.z)
	pi.1 = 1-pi.z-pi.0
	
	if(verbose) message("starting DEup.mean...")
	u = find.mean.up(.99, null.mean = fdrmodel$fp0[3,1], null.sd = fdrmodel$fp0[3,2], null.prop = pi.0.tmp, vals = nzt)
	DEup.mean = u$m
	if(verbose) message("starting DEdown.mean...")
	d = find.mean.down(.01, null.mean = fdrmodel$fp0[3,1], null.sd = fdrmodel$fp0[3,2], null.prop = pi.0.tmp, vals = nzt)
	DEdown.mean = d$m
	DEup.sd = find.sd(prev.p = u$p, found.mean = DEup.mean, null.mean = fdrmodel$fp0[3,1], null.sd = fdrmodel$fp0[3,2], null.prop = pi.0.tmp, vals = nzt, up=TRUE)
	DEdown.sd = find.sd(prev.p = d$p, found.mean = DEdown.mean, null.mean = fdrmodel$fp0[3,1], null.sd = fdrmodel$fp0[3,2], null.prop = pi.0.tmp, vals = nzt, up=FALSE)

	if(plots){
		zzaxis <- seq(-10,10,length=2000)
		hist(tstats,breaks=100,col="gray",freq=F,ylim=c(0,0.4),main="fitted distributions",xlim=c(-5,5),xlab="t statistics")
		lines(zzaxis,pi.0*dnorm(zzaxis,fdrmodel$fp0[3,1],fdrmodel$fp0[3,2]),lwd=3)
		lines(zzaxis,(pi.1)/2*dnorm(zzaxis,DEup.mean,DEup.sd),lwd=3,col="red")
		lines(zzaxis,(pi.1)/2*dnorm(zzaxis,DEdown.mean,DEdown.sd),lwd=3,col="green")
		legend("topright",lty=c(1,1,1),lwd=c(3,3,3),col=c("black","red","green"),c("null","DE up","DE down"))
		if(!is.null(plotfile)) dev.off()
	}
	
	params = list(mean=c(0,fdrmodel$fp0[3,1],DEup.mean,DEdown.mean),sd=c(0.0000001,fdrmodel$fp0[3,2],DEup.sd,DEdown.sd))
	out = list(params=params, stateprobs=c(pi.z,pi.0,pi.1/2,pi.1/2),pi.z=pi.z, pi.0=pi.0, pi.1=pi.1)
	return(out)
}

#test:
load("~/Documents/hopkins/research/_readlets/_zeiger-readlets/st22.rda")
#getParams(st,plots=F,plotfile="~/Documents/hopkins/research/_ebHMM-Rfunctions/testplot.pdf")
getParams(st,plots=T)

# test on different data:
#load("~/Documents/hopkins/research/_sarven-brainSamples/st/st22.rda")
#getParams(st,plots=T,plotfile="~/Documents/hopkins/research/_ebHMM-Rfunctions/testplot-sarven.pdf")




