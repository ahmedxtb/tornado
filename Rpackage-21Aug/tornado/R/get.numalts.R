get.numalts <-
function(pctil,null.mean,null.sd,null.prop,vals,up = TRUE){
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
