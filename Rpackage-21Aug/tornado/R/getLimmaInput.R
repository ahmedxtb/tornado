getLimmaInput <-
function(dbfile, tablename, group, chunksize = 100000, adjustvars = NULL, colsubset = NULL){
	require(limma)
	require(multicore)
	require(Genominator)
	
	# get ready to read from database:
	tab = ExpData(dbfile, tablename)
	if(!is.null(colsubset)) tab = tab[,colsubset] #not recommended, as this loads matrix into memory.
	pos = tab[,1]$pos
	N = length(pos)
	lastloop = trunc(N/chunksize)

	# create model matrix:
	ncol = length(tab[1,])
	colmeds = NULL
	for(i in 2:ncol){
		eval(parse(text=paste("med = median(tab[,",i,"]$",names(tab[,i]),")",sep="")))
		colmeds[i-1] = med
	} #get median of each column to use as adjustment variable
	if(!is.null(adjustvars)){
		string1 = ""
		for(i in 1:dim(adjustvars)[2]){
			eval(parse(text=paste("av",i," <- adjustvars[,",i,"]",sep="")))
			string1 = paste(string1, paste("av",i,sep=""),sep="+")
		}
		eval(parse(text=paste("x = model.matrix(~group+colmeds",string1,")",sep="")))
	}else{x = model.matrix(~group+colmeds)}	
	
	# define modeling function to apply to each chunk:
	lmFit.apply = function(i){
  		if(i!=lastloop) mymat <- tab[(chunksize*i+1):(chunksize*(i+1)),-1] #-1 removes pos
  		else mymat <- tab[(chunksize*i+1):N,-1] #-1 removes pos
  		mymat <- log2(mymat+0.5)
  		Amean <- rowMeans(mymat) 
  		fit <- lmFit(mymat,x)
  		return(list(fit=fit, Amean=Amean))
  		}
	
	# fit a model to each row (chunk) of database:
	if(interactive()) {
		warning("Cannot use mclapply in interactive session; reverting to single-core lapply")
		lmFit.output = lapply(0:lastloop, lmFit.apply)
		}
	if(!interactive()) {
		warning("mclapply functionality not yet implemented: fitting models with one core.")
		lmFit.output = lapply(0:lastloop, lmFit.apply)
	}
	# gather results from all the models and return them (this part isn't tested):
	coef = stdev = sma = dfres = am = NULL
	for(i in 1:length(lmFit.output)){
  		coef = append(coef,lmFit.output[[i]]$fit$coefficients[,2])
  		stdev = append(stdev, lmFit.output[[i]]$fit$stdev.unscaled[,2])
  		sma = append(sma, lmFit.output[[i]]$fit$sigma)
  		dfres = append(dfres, lmFit.output[[i]]$fit$df.residual)
  		am = append(am, lmFit.output[[i]]$Amean)
	}
	return(list(ebobject = list(coefficients = as.numeric(coef), stdev.unscaled = as.numeric(stdev), sigma = sma, df.residual = dfres, Amean = am), pos = pos))
	
}
