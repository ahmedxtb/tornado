find.sd <-
function(prev.p, found.mean, null.mean, null.sd, null.prop, vals, up=T){
	if(up) new.percentile = (1+prev.p)/2 # we increase the percentile we're going to look at
	if(!up) new.percentile = prev.p/2 #decrease it.
	cutoff.val = qnorm(new.percentile, null.mean, null.sd) # x-value that is that percentile of the theoretical null distribution
	num.nulls = null.prop*length(vals)
	coef1 = ifelse(up,1-new.percentile,new.percentile)
	num.nulls.above = coef1*num.nulls
	num.total.above = ifelse(up,sum(vals>cutoff.val),sum(vals<cutoff.val))
	num.alts.above = num.total.above - num.nulls.above
	num.alts = (1-null.prop)*0.5*length(vals) # this is the number of DEup values
	if(up) alt.percentile = 1-(num.alts.above/num.alts)
	if(!up) alt.percentile = num.alts.above/num.alts
	zstat = qnorm(alt.percentile)
	if(is.na(zstat)){
		warning("Numerical standard deviation estimation failed.  Defaulting to sd of estimated null distribution.")
		gp = getParams.failsafe(null.mean, null.sd)
		if(up) return(gp$DEup.sd)
		if(!up) return(gp$DEdown.sd)
	}
	sigma = (cutoff.val - found.mean)/zstat
    if(num.alts.above<=0 | sigma<=0) {
		warning("Numerical standard deviation estimation failed. Defaulting to sd of estimated null distribution.")
		gp = getParams.failsafe(null.mean, null.sd)
		if(up) return(gp$DEup.sd)
		if(!up) return(gp$DEdown.sd)
	}
	return(sigma)
}
