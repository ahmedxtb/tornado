getParams.failsafe <-
function(null.mean, null.sd){
		DEup.mean = qnorm(.95, mean=null.mean, sd=null.sd)
		DEup.sd = DEdown.sd = null.sd
		DEdown.mean = qnorm(.05, mean=null.mean, sd=null.sd)
		return(list(DEup.mean = DEup.mean, DEdown.mean = DEdown.mean, DEup.sd = DEup.sd, DEdown.sd = DEdown.sd))
}
