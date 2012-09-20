getSvs <- function(dbfile, tablename, group, chunksize = 100000,colsubset = NULL){
  require(limma)
  require(multicore)
  require(Genominator)
  require(sva)
  require(genefilter)
	
  tab = ExpData(dbfile, tablename)
  pos = as.vector(tab[,1])
  if(!is.null(colsubset)) tab = tab[,c(1,colsubset)] #not recommended, as this loads matrix into memory.
  print("subset")
  
  N = dim(pos)[1]
  lastloop = trunc(N/chunksize)
  print(lastloop)
  
  vv.apply = function(i){
    if(i!=lastloop) mymat <- tab[(chunksize*i+1):(chunksize*(i+1)),-1] #-1 removes pos
    else mymat <- tab[(chunksize*i+1):N,-1] #-1 removes pos
    mymat <- log2(as.matrix(mymat)+0.5)
    rowV <- rowVars(mymat)
    return(list(rowV=rowV))
  }

  vv.output = lapply(0:lastloop,vv.apply)

  # Collect all the variances
  vv = NULL
  for(i in 1:length(vv.output)){
    vv = append(vv,vv.output[[i]]$rowV)  
  }
  print(length(vv))
  ind = which(rank(-vv) <= 10000)
  mod = model.matrix(~as.factor(group))
  mod0 = cbind(mod[,1])
  print(as.matrix(tab[ind,-1])[1,])
  svaObj = sva(as.matrix(tab[ind,-1]),mod=mod,mod0=mod0)
  return(svaObj)
}




