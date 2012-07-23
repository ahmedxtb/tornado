## updated 7/23/12 AF
## database file (containing filtered coverage) --> precursors to eBayes/getting moderated t-statistics.
## i.e. in pipeline: this goes after makeDb.R, before getTstats.R

getLimmaInput <- function(dbfile, tablename, adjustvars = NULL, group, chunksize = 100000, colsubset = NULL){
	#dbfile and tablename are from table created with makeDb()
	#adjustvars is either left null or is a matrix (nxp) of adjustment variables (SVA?)
	#group is a length n 0/1 vector grouping the samples
	
	require(Genominator)
	require(limma)
	tab = ExpData(dbfile, tablename)
	pos = tab[,1]$pos #ASSUMES FIRST COLUMN IS CALLED POS!!!!!
	N = length(pos)
	
	if(!is.null(colsubset)) tab = tab[,colsubset]
	
	# set up model matrix (we are regressing count on group plus optional adjustment variables)
	if(!is.null(adjustvars)){
		string1 = ""
		for(i in 1:dim(adjustvars)[2]){
 		  eval(parse(text=paste("av",i," <- adjustvars[,",i,"]",sep="")))
		  string1 = paste(string1,paste("av",i,sep=""),sep="+")
		}
		eval(parse(text=paste("x = model.matrix(~group",stringarg,")",sep="")))
	}else{x = model.matrix(~group)}
	
	# set up empty vectors we need to fill
	listcoefficients = matrix(NA,nrow=N,ncol=dim(x)[2])
	liststdev.unscaled = matrix(NA,nrow=N,ncol=dim(x)[2])
	listsigma = listAmean = rep(NA,N)
	listdf.residual <- rep(length(tab[1,])-dim(x)[2],N) #no need to get this every time.

	# loop through database to create these vectors - these are input to eBayes
	lastloop = trunc(N/chunksize)
	for(i in 0:lastloop){
		if(i!=lastloop) mymat <- tab[(chunksize*i+1):(chunksize*(i+1)),]
		else mymat <- tab[(chunksize*i+1):N,]
		mymat <- log2(mymat+0.5)
	    mymat <- sweep(mymat,2,Biobase::rowMedians(t(mymat)))
		fit <- lmFit(mymat,x)
		startind <- which(is.na(listsigma)==TRUE)[1]
		endind <- startind+dim(mymat)[1]-1
		inds <- c(startind:endind)
		listcoefficients[inds,] <- fit$coefficients
		liststdev.unscaled[inds,] <- fit$stdev.unscaled
		listsigma[inds] <- fit$sigma
		listAmean[inds] <- rowMeans(mymat)
		rm(mymat,fit);gc();gc()
		}
	
	return(list(ebobject = list(coefficients = listcoefficients[,2], stdev.unscaled = liststdev.unscaled[,2], sigma = listsigma, df.residual = listdf.residual, Amean = listAmean), pos = pos))
	
} #end function

# recommendation: pass ebobject to getTstats(), save pos as rda file to use later.

