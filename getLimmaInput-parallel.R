# parallelized version of getLimmaInput - runs faster
# AF 8/15/12
# NOT YET TESTED AS FUNCTION (though the pieces inside have been run)


## getLimmaInput.p()
## arguments:
## --dbfile: string giving location of .db file (usually created with makeDb()) containing the data
## --tablename: name of table in the .db file (usually created with makeDb())
## (note: dbfile and tablename arguments in makeDb() and this function should match)
## --adjustvars: if desired, an nxp matrix of adjustment variables (confounders, SVA output, etc.)
## --group: a length n 0/1 vector grouping the samples (currently only handles 2 groups)
## --colsubset: if desired, column indices of the input file of the samples you wish to include in analysis. Should NOT include 1 (pos).
## --chunksize: sections of the database are analyzed separately - how big (how many rows) should those pieces be? 
## return:
## a list containing elements $ebobject, which is a list mimicking the output of running lmFit on the whole dataset,
## and $pos, giving the row indices of all the positions on the chromosome passing the filtering criterion (i.e. bases with nontrivial coverage)
## usage recommendataion: run getLimmaInput.p(), pass $ebobject to getTstats(), save $pos as rda file in case it is needed later

getLimmaInput.p <- function(dbfile, tablename, adjustvars, group, colsubset, chunksize = 100000){
	require(limma)
	require(multicore)
	require(Genominator)
	
	# create model matrix:
	if(!is.null(adjustvars)){
		string1 = ""
		for(i in 1:dim(adjustvars)[2]){
			eval(parse(text=paste("av",i," <- adjustvars[,",i,"]",sep="")))
			string1 = paste(string1, paste("av",i,sep=""),sep="+")
		}
		eval(parse(text=paste("x = model.matrix(~group",string1,")",sep="")))
	}else{x = model.matrix(~group)}	
	
	# get ready to read from database:
	tab = ExpData(dbfile, tablename)
	pos = tab[,1]$pos
	N = length(pos)
	lastloop = trunc(N/chunksize)
	
	# define modeling function to apply to each chunk:
	lmFit.apply = function(i){
  		if(i!=lastloop) mymat <- tab[(chunksize*i+1):(chunksize*(i+1)),-1] #-1 removes pos
  		else mymat <- tab[(chunksize*i+1):N,-1] #-1 removes pos
  		mymat <- log2(mymat+0.5)
  		Amean <- rowMeans(mymat) # is this the mean we want? or do we want the mean of the centered values? or a different adjustment?
  		mymat <- sweep(mymat,2,Biobase::rowMedians(t(mymat)))
  		fit <- lmFit(mymat,x)
  		return(list(fit=fit, Amean=Amean))
  		}
	
	# fit a model to each row (chunk) of database:
	if(interactive()) {
		warning("Cannot use mclapply in interactive session; reverting to single-core lapply")
		lmFit.output = lapply(1:lastloop, lmFit.apply)
		}
	if(!interactive()) lmFit.output = mclapply(1:lastloop, lmFit.apply)
	
	# gather results from all the models and return them (this part isn't tested):
	coef = stdev = sma = dfres = am = NULL
	for(i in 1:length(lmFit.output)){
  		coef = append(coef,lmFit.output[[i]]$fit$coefficients[,2])
  		stdev = append(stdev, lmFit.output[[i]]$fit$stdev.unscaled[,2])
  		sma = append(sma, lmFit.output[[i]]$fit$sigma)
  		dfres = append(dfres, lmFit.output[[i]]$fit$df.residual)
  		am = append(am, lmFit.output[[i]]$Amean)
	}
	return(list(ebobject = list(coefficients = coef, stdev.unscaled = stdev, sigma = sma, df.residual = dfres, Amean = am), pos = pos))
	
} #end function
