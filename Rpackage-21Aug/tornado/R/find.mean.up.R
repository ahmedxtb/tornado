find.mean.up <-
function(init.value, null.mean, null.sd, null.prop, vals){
	if(init.value<=0.5) stop("finding DE up quantile - init.value should be >0.5")
	x = init.value
	target = round(0.5*0.5*(1-null.prop)*length(vals))  #half the expected # of DEup values.
	history = pseq = numseq = qseq = NULL
	counter = 0
	while(TRUE){
		iter = try(get.numalts(x,null.mean,null.sd,null.prop,vals,up=T),silent=T)
		if(class(iter)=="try-error") {
			warning("Numerical estimation of DE-up mean failed.  Defaulting to mean = 95th percentile of estimated null distribution, sd = sd of estimated null distribution.")
			gp = getParams.failsafe(null.mean, null.sd)
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
				warning("Numerical estimation of DE-up mean failed.  Defaulting to mean = 95th percentile of estimated null distribution, sd = sd of estimated null distribution.")
				gp = getParams.failsafe(null.mean, null.sd)
				return(list(m=gp$DEup.mean, s=gp$DEup.sd))
			}			
			if(history[counter-1]==-1 & numseq[counter]<=numseq[counter-1]){
				warning("Numerical estimation of DE-up mean failed.  Defaulting to mean = 95th percentile of estimated null distribution, sd = sd of estimated null distribution.")
				gp = getParams.failsafe(null.mean, null.sd)
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
