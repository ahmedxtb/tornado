# helper functions for getParams.R
# 7/18/12

last <- function(x) return(tail(x, n=1))

get.numalts <- function(pctil,null.mean,null.sd,null.prop,vals,up = TRUE){
	cutoff.val = qnorm(pctil,null.mean,null.sd)
	num.nulls = null.prop*length(vals)
	# I use "above" to mean "above or below". very loose...
	coef1 = ifelse(up,1-pctil,pctil)
	num.nulls.above = coef1*num.nulls
	num.total.above = ifelse(up,sum(vals>cutoff.val),sum(vals<cutoff.val))
	num.alts.above = num.total.above - num.nulls.above
	if(num.alts.above<=0) stop("negative num.alts.above. rethink algorithm. sad day.")
	return(list(num=round(num.alts.above),val=cutoff.val))
}

find.mean.up <- function(init.value, null.mean, null.sd, null.prop, vals, verbose=F){
	if(init.value<=0.5) stop("finding DE up quantile - init.value should be >0.5")
	x = init.value
	target = round(0.5*0.5*(1-null.prop)*length(vals))  #half the expected # of DEup values.
	history = pseq = NULL
	counter = 0
	while(TRUE){
		iter = get.numalts(x,null.mean,null.sd,null.prop,vals,up=T)
		if(iter$num == target) return(list(m=iter$val,p = pseq[counter]))
		if(counter>1){
			if(abs(pseq[counter]-pseq[counter-1])<1e-8) return(list(m=iter$val,p=pseq[counter]))
		}
		counter = counter+1
		if(verbose) print(counter)
		pseq[counter] = x
		if(verbose) print(x)
		
		if(iter$num > target){
			history[counter] = 1
			if(counter==1){ x = (x+1)/2; next }
		}
		if(iter$num < target){
			history[counter] = -1
			if(counter==1){ x = (x+0.5)/2; next }
		}
		
		if(history[counter-1] != history[counter]) x = (x+pseq[counter-1])/2
		if(history[counter-1] == history[counter]){
				prev = last(which(history != history[counter]))
				if(length(prev)==0 & iter$num > target){x = (x+1)/2; next}
				if(length(prev)==0 & iter$num < target){x = (x+0.5)/2; next}
				x = (x+pseq[prev])/2
			}
		}
}

find.mean.down <- function(init.value, null.mean, null.sd, null.prop, vals, verbose=F){
	if(init.value>=0.5) stop("finding DE down quantile - init.value should be <0.5")
	x = init.value
	target = round(0.5*0.5*(1-null.prop)*length(vals))  #half the expected # of DEdown values.
	history = pseq = NULL
	counter = 0
	while(TRUE){
		iter = get.numalts(x,null.mean,null.sd,null.prop,vals,up=F)
		if(iter$num == target) return(list(m=iter$val,p = pseq[counter]))
		if(counter>1){
			if(abs(pseq[counter]-pseq[counter-1])<1e-8) return(list(m=iter$val,p=pseq[counter]))
		}
		counter = counter+1
		if(verbose) print(counter)
		pseq[counter] = x
		if(verbose) print(x)
		
		if(iter$num > target){
			history[counter] = 1
			if(counter==1){ x = x/2; next }
		}
		if(iter$num < target){
			history[counter] = -1
			if(counter==1){ x = (x+0.5)/2; next }
		}
		
		if(history[counter-1] != history[counter]) x = (x+pseq[counter-1])/2
		if(history[counter-1] == history[counter]){
				prev = last(which(history != history[counter]))
				if(length(prev)==0 & iter$num > target){x = x/2; next}
				if(length(prev)==0 & iter$num < target){x = (x+0.5)/2; next}
				x = (x+pseq[prev])/2
			}
		}
}

find.mean <- function(init.value, null.mean, null.sd, null.prop, vals, up = TRUE){
	if(up) k = find.mean.up(init.value, null.mean, null.sd, null.prop, vals)
	if(!up) k = find.mean.down(init.value, null.mean, null.sd, null.prop, vals)
	return(k)
}


find.sd <- function(prev.p, found.mean, null.mean, null.sd, null.prop, vals, up=T){
	if(up) new.percentile = (1+prev.p)/2 # we increase the percentile we're going to look at
	if(!up) new.percentile = prev.p/2 #decrease it.
	cutoff.val = qnorm(new.percentile, null.mean, null.sd) # x-value that is that percentile of the theoretical null distribution
	num.nulls = null.prop*length(vals)
	coef1 = ifelse(up,1-new.percentile,new.percentile)
	num.nulls.above = coef1*num.nulls
	num.total.above = ifelse(up,sum(vals>cutoff.val),sum(vals<cutoff.val))
	num.alts.above = num.total.above - num.nulls.above
	if(num.alts.above<=0) stop("yucky! num.alts.above<0!")
	num.alts = (1-null.prop)*0.5*length(vals) # this is the number of DEup values
	if(up) alt.percentile = 1-(num.alts.above/num.alts)
	if(!up) alt.percentile = num.alts.above/num.alts
	zstat = qnorm(alt.percentile)
	sigma = (cutoff.val - found.mean)/zstat
	return(sigma)
}




