find.mean <-
function(init.value, null.mean, null.sd, null.prop, vals, up = TRUE){
	if(up) k = find.mean.up(init.value, null.mean, null.sd, null.prop, vals)
	if(!up) k = find.mean.down(init.value, null.mean, null.sd, null.prop, vals)
	return(k)
}
