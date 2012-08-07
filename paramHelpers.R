# helper functions for getParams.R
# 7/18/12
# updated 8/7/12 to include the failsafe function

last <- function(x) return(tail(x, n=1))

getParams.failsafe <- function(){
	# this shouldn't need any arguments, since it's used entirely within the other functions.
	# used if the other method fails.
			if(!exists("fdrmodel")) stop("fdrmodel must be defined")
			bincounts = fdrmodel$yt
			binlocs = fdrmodel$x
			downinds = which(binlocs<0)
			upinds = which(binlocs>=0)
			oecx = binlocs[upinds]
			oecy = bincounts[upinds]
			uecx = binlocs[downinds]
			uecy = bincounts[downinds]
			mystats <- mystats2 <- list()
			for(i in 1:length(oecx)){
                          if(oecy[i]==0) mystats[[i]] <- NA
                          if(oecy[i]!=0) mystats[[i]] <- rep(oecx[i],round(oecy[i]))
			}
			mystats <- unlist(mystats)
			DEup.sd <- sd(mystats,na.rm=TRUE)
			DEup.mean <- mean(mystats,na.rm=TRUE)
			for(i in 1:length(uecx)){
                          if(uecy[i]==0) mystats2[[i]] <- NA
                          if(uecy[i]!=0) mystats2[[i]] <- rep(uecx[i],round(uecy[i]))
			}
			mystats2 <- unlist(mystats2)
			DEdown.sd <- sd(mystats2,na.rm=TRUE)
			DEdown.mean <- mean(mystats2,na.rm=TRUE)
			return(list(DEup.mean = DEup.mean, DEdown.mean = DEdown.mean, DEup.sd = DEup.sd, DEdown.sd = DEdown.sd))
}

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

find.mean.up <- function(init.value, null.mean, null.sd, null.prop, vals){
	if(init.value<=0.5) stop("finding DE up quantile - init.value should be >0.5")
	x = init.value
	target = round(0.5*0.5*(1-null.prop)*length(vals))  #half the expected # of DEup values.
	history = pseq = numseq = qseq = NULL
	counter = 0
	while(TRUE){
		iter = try(get.numalts(x,null.mean,null.sd,null.prop,vals,up=T),silent=T)
		if(class(iter)=="try-error") {
			warning("Numerical estimation of DE-up mean failed.  Using locfdr estimates.")
			gp = getParams.failsafe()
			return(list(m=gp$DEup.mean, s=gp$DEup.sd))
		}
		if(iter$num == target) return(list(m=iter$val,p = pseq[counter]))
		if(counter>1){
			if(abs(pseq[counter]-pseq[counter-1])<1e-8) return(list(m=iter$val,p=pseq[counter]))
		}
		counter = counter+1
		numseq[counter] = iter$num
		qseq[counter] = iter$val
		pseq[counter] = x
		
		if(counter>1){
			if(history[counter-1]==1 & numseq[counter]>=numseq[counter-1]){
				warning("numerical estimation of DE-up mean failed.  Using locfdr estimates.")
				gp = getParams.failsafe()
				return(list(m=gp$DEup.mean, s=gp$DEup.sd))
			}
			if(history[counter-1]==-1 & numseq[counter]<=numseq[counter-1]){
				warning("numerical estimation of DE-up mean failed.  Using locfdr estimates.")
				gp = getParams.failsafe()
				return(list(m=gp$DEup.mean, s=gp$DEup.sd))
			}
		}
		
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

find.mean.down <- function(init.value, null.mean, null.sd, null.prop, vals){
	if(init.value>=0.5) stop("finding DE down quantile - init.value should be <0.5")
	x = init.value
	target = round(0.5*0.5*(1-null.prop)*length(vals))  #half the expected # of DEdown values.
	history = pseq = numseq = qseq = NULL
	counter = 0
	while(TRUE){
		iter = try(get.numalts(x,null.mean,null.sd,null.prop,vals,up=F),silent=T)
		if(class(iter)=="try-error") {
			warning("numerical estimation of DE-down mean failed. Using locfdr estimates.")
			gp = getParams.failsafe()
			return(list(m=gp$DEdown.mean,s=gp$DEdown.sd))
		}
		if(iter$num == target) return(list(m=iter$val,p = pseq[counter]))
		if(counter>1){
			if(abs(pseq[counter]-pseq[counter-1])<1e-8) return(list(m=iter$val,p=pseq[counter]))
		}
		counter = counter+1
		numseq[counter] = iter$num
		qseq[counter] = iter$val
		pseq[counter] = x
		
		if(counter>1){
			if(history[counter-1]==1 & numseq[counter]>=numseq[counter-1]){
				warning("numerical estimation of DE-down mean failed.  Using locfdr estimates.")
				gp = getParams.failsafe()
				return(list(m=gp$DEdown.mean, s=gp$DEdown.sd))
			}
			if(history[counter-1]==-1 & numseq[counter]<=numseq[counter-1]){
				warning("numerical estimation of DE-down mean failed.  Using locfdr estimates.")
				gp = getParams.failsafe()
				return(list(m=gp$DEdown.mean, s=gp$DEdown.sd))
			}
		}
		
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
	num.alts = (1-null.prop)*0.5*length(vals) # this is the number of DEup values
	if(up) alt.percentile = 1-(num.alts.above/num.alts)
	if(!up) alt.percentile = num.alts.above/num.alts
	zstat = qnorm(alt.percentile)
	if(is.na(zstat)){
		warning("Numerical standard deviation estimation failed.  Using locfdr estimates.")
		gp = getParams.failsafe()
		if(up) return(gp$DEup.sd)
		if(!up) return(gp$DEdown.sd)
	}
	sigma = (cutoff.val - found.mean)/zstat
    if(num.alts.above<=0 | sigma<=0) {
		warning("Numerical standard deviation estimation failed. Using locfdr estimates.")
		gp = getParams.failsafe()
		if(up) return(gp$DEup.sd)
		if(!up) return(gp$DEdown.sd)
	}
	return(sigma)
}







