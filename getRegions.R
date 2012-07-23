# function going from estimated distribution params --> list of regions
# alyssa frazee 6/6/12

getRegions <- function(method, chromosome, pos, tstats, transprobs = c(0.999, 1e-12), stateprobs = NULL, 
	params = NULL, K = 25, tcut = 2, includet=F, includefchange=F, fchange=NULL)
{
  	if(method!="HMM" & method!="smoothcut" & method!="CBS") stop("Invalid method. Choices are HMM, smoothcut, or CBS")
  	if(sum(sort(pos)!=pos)>0) stop("pos (and probably t-statistics) improperly sorted")
    if (method == "HMM"){
		require(HiddenMarkov)
		if(is.null(params)) stop("HMM method requires parameter estimation - please provide theoretical values or estimate with getParams")
		stayprob = transprobs[1]
    	EtoDE = transprobs[2]
    	EtoZero = 1-stayprob-2*EtoDE
    	transmat = matrix(c(stayprob,(1-stayprob)/3,(1-stayprob)/3,(1-stayprob)/3,				 							         							EtoZero,stayprob,EtoDE,EtoDE,
					     	EtoZero,EtoDE,stayprob,EtoDE,
					     	EtoZero,EtoDE,EtoDE,stayprob),nrow=4,ncol=4,byrow=T)
    	if(is.null(stateprobs)) stop("HMM method requires initial values for each state - please provide theoretical values or estimate with getParams")
    	if(sum(names(params)==c("mean","sd")) < length(params)) stop("params argument formatted incorrectly; please see getParams")
    
   		hmmodel = dthmm(x=tstats, Pi=transmat, delta=stateprobs, distn="norm",pm=params) 
    	state.predictions = Viterbi(hmmodel)
		} #end HMM
		
	if (method == "smoothcut"){	
		smootht = runmed(tstats,k=K,endrule="constant")
		state.predictions = rep(NA,length(smootht))
		state.predictions[which(smootht==0)] = 1
		state.predictions[which(smootht>0 & smootht<tcut)] = 2
		state.predictions[which(smootht<0 & smootht>-tcut)] = 2
		state.predictions[which(smootht>=tcut)] = 3
		state.predictions[which(smootht<=-tcut)] = 4
		} # end smoothcut
		
	if(method == "CBS"){	
		if(is.null(params)) stop("CBS method requires parameter estimation - please provide theoretical values or estimate with getParams")
		require(DNAcopy)
		cna.object = CNA(tstats,chrom=rep(chromosome,length(tstats)),maploc=pos,data.type="logratio",sampleid="all",presorted=TRUE)
		cna.smoothed = smooth.CNA(cna.object)
		segmentation = segment(cna.smoothed,p.method="hybrid",verbose=0) # slow.

		statecalls = segmentation$output$seg.mean

		zprobs = 1-pnorm(abs((statecalls-params$mean[1])/params$sd[1]))
    	nullprobs = 1-pnorm(abs((statecalls-params$mean[2])/params$sd[2]))
		upprobs = 1-pnorm(abs((statecalls-params$mean[3])/params$sd[3]))
		downprobs = 1-pnorm(abs((statecalls-params$mean[4])/params$sd[4]))
		state.predictions.segments = max.col(cbind(zprobs,nullprobs,upprobs,downprobs))
		statetable = list(lengths=segmentation$output$num.mark, values=state.predictions.segments)
		state.predictions = inverse.rle(statetable)
		} # end CBS
		
    states.norle = data.frame(chr=rep(chromosome,length(pos)),pos,state=state.predictions)
    
    pos.full = c(pos[1]:pos[length(pos)])
    state.preds.full = rep(1,length(pos.full)) 
    state.preds.full[pos-pos[1]+1] = state.predictions
    full.rle = rle(state.preds.full)
    END = pos.full[cumsum(full.rle$lengths)]
    START = END-full.rle$lengths+1
    states = data.frame(chr=rep(chromosome,length(START)),start=START,end=END,state=full.rle$values,length=full.rle$lengths)
    if(includet & !includefchange){
    	full.t = rep(0,max(pos))
    	full.t[pos] = tstats
    	tmean = NULL
    	for(i in 1:dim(states)[1]){
    		tmean <- mean(full.t[states$start[i]:states$end[i]])
    	}
    	states = data.frame(states,mean.t=tmean)
    }
    if(!includet & includefchange){
    	if(is.null(fchange)) stop("must provide fold changes (log2 scale) if you want to include them in your output (output fold change is not on log scale)")
    	full.fchange = rep(1,max(pos))
    	full.fchange[pos] = 2^fchange
    	fcmean = NULL
    	for(i in 1:dim(states)[1]){
    		fcmean <- mean(full.fchange[states$start[i]:states$end[i]])
    	}
    	states = data.frame(states,mean.fold.change=fcmean)
    }
	if(includet & includefchange){
    	full.t = rep(0,max(pos))
    	full.t[pos] = tstats
    	full.fchange = rep(1,max(pos))
    	full.fchange[pos] = 2^fchange
    	tmean = fcmean = NULL
    	for(i in 1:dim(states)[1]){
       		fcmean <- mean(full.fchange[states$start[i]:states$end[i]])
    		tmean <- mean(full.t[states$start[i]:states$end[i]])
    	}
    	states = data.frame(states,mean.t=tmean,mean.fold.change=fcmean)
	}
    out = list(states.norle=states.norle, states=states)
    return(out)

} #end function






