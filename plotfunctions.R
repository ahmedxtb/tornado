## plotting functions for HMM (or whatever) results
## ACF June 8 2012
## includes plotRegion, plotExon, and plotGene

# currently, not a lot of plotting arguments are available here - add flexibility later?
# only tested with already-loaded matrix (haven't tried reading from file or database yet) 7.23.12

## plotRegion():
## arguments:
## --getRegionObject: output from getRegions (list contatining $states and $states.norle)
## --ind: index of the region in getRegionObject$states that you're interested in plotting
## --tstats: vector of t-statistics used in getRegions
## --pos: vector of positions corresponding to tstats
## --annotation: data frame containing exon annotation to use (see getAnnotation)
## --counts: raw data used in getting the t-statistics. this can be either a string indicating the .db file created with makeDb, 
## --counts(cont): ...a string indicating the location of a text file containing coverage (probably not a good option, since the raw matrix probably won't have the zero-rows removed)
## --counts(cont): ...or an already-loaded matrix containing the raw data.
## --counts(cont): NOTE that count must have the same number of rows as the # of elements in tstats and pos, and the rows must correspond to pos.
## --group: a vector containing the group labels for the columns of counts. Only 2 groups are permitted at this time.
## --bppad: the number of bases to plot outside of the designated region (default 50).  Essentially, use this to "zoom" in (decrease bppad) or out (increase bppad) on the region.
## --axpad: how much wider (in bases) you'd like the x-axis to be, compared to the plotted area
## --prettyskips: if TRUE, plot counts/states/t-statistics contiguously, even if there are zero entries between them (i.e., even though pos may not indicate that contiguous postions are being plotted).
## --prettyskips(cont): note that in general, when plotting just one region, this is not an issue as regions tend to be contiguous. Will only affect areas outside the region, i.e., has larger impact if bppad is large.
## --skiplines: if TRUE, add a light vertical line to the plot indicating an eliminated "zero" region
## --countsheader: if TRUE, the counts matrix contains a header row. Not usually the case if counts is a database or already-loaded matrix.
## --countssep: (defaults to tab) - if reading counts from a text file, the separator used in that file.
## --tabname: if counts is a database file, the name of the table that was dumped into that database.
## --plotfile: optional string containing a file name you'd like to put the plot into (if NULL, plot appears interactively). Should be a .jpg file at this time.
## --width and height: only used with plotfile - dimensions (in pixels) of resulting jpg.  defaults to 900x750
## --plottitle: optional main title to use on your plot.  Defaults to chromosome: start-end (referring to plotted REGION)
## --chromosome: the chromosome corresponding to the region you're plotting, IN THE SAME FORMAT AS IS INCLUDED IN THE SUPPLIED annotation ARGUMENT
## --legendloc: string indicating one of "topright","bottomright","topleft",or "bottomleft" indicating where the legend (indicating group label on raw count plot) should be located. Defaults to "bottomleft."
## return:
## a lovely plot of a specified region. (either interactively or in a jpeg file)

plotRegion <- function(getRegionObject, ind, tstats, pos, annotation, counts, group, bppad = 50, axpad = 50, prettyskips = T, skiplines = T, countsheader=F, countssep="\t", tabname = NULL, plotfile = NULL, width = 900, height = 750, plottitle = NULL, chromosome,legendloc="bottomleft")
{
	if(!all(c("states.norle","states") %in% names(getRegionObject))) stop("getRegionObject must be a list with elements \"states.norle\" and \"states\" - usually the return of getRegions()")
	if(!all(c("chr","start","end") %in% names(getRegionObject$states))) stop("getRegionObject \"states\" component must contain chr, start, and end")
	chr = getRegionObject$states$chr[ind]
	Start = getRegionObject$states$start[ind]
	End = getRegionObject$states$end[ind]
	if(!all(c("chr","start","end")%in%colnames(annotation))) stop("annotation must contain columns named chr, start, and end")
	if(sum(sort(pos)!=pos)>0) stop("pos improperly sorted")
	if(!prettyskips) stop("I haven't implemented this without prettyskips yet - sorry!")
	
	if(is.character(counts)){
		if(substr(counts,nchar(counts)-2,nchar(counts)==".db")){
			require(Genominator)
			if(is.null(tabname)) stop("reading matrix from database requires table name")
			counts <- ExpData(dbFilename = counts, tablename = tabname)
		}
		if(substr(counts,nchar(counts)-2,nchar(counts)!=".db")){
		warning("Reading count matrix from file. May take several minutes and/or exceed memory limits.")
		counts <- read.table(counts,header=countsheader,sep=countssep)
		}
	}
	
	# other choice: counts can be an already-loaded matrix
	
	lowerbound = Start-bppad
	upperbound = End+bppad
	xaxinds = which(pos>=lowerbound & pos<=upperbound)
	plotlb = xaxinds[1]-axpad
	plotub = xaxinds[length(xaxinds)]+axpad
	
	if(!is.null(plotfile)) jpeg(file=plotfile,width=width,height=height)
	
	par(mar=c(1,3.2,2,1),mgp=c(1.5,0.5,0),mfrow=c(3,1),cex.lab=1.5,omi=c(0.4,0,0,0))
	
	## panel A: plots of raw data
	groupsplit = split(1:length(group),group)
	if(length(groupsplit)!=2) stop("need exactly 2 comparison groups")
	g1 = groupsplit[[1]]
	g2 = groupsplit[[2]]
	firstcolor = ifelse(1 %in% groupsplit[[1]],"blue","red")
	plot(xaxinds,log2(counts[,1]+0.5)[xaxinds],type="l",lwd=0.5,col=firstcolor,ylim=c(-2,9),xaxt="n",
		xlim=c(plotlb,plotub),xlab="",ylab="log2 count")
	for(i in g1){
  		lines(xaxinds,log2(counts[xaxinds,i]+0.5),type="l",lwd=0.5,col="blue")
	}
	for(i in g2){
  		lines(xaxinds,log2(counts[xaxinds,i]+0.5),type="l",lwd=0.5,col="red")
	}
	mean1 = rowMeans(log2(counts[xaxinds,g1]+0.5))
	mean2 = rowMeans(log2(counts[xaxinds,g2]+0.5))
	lines(xaxinds,mean1,col="blue",lwd=3)
	lines(xaxinds,mean2,col="red",lwd=3)
	if(is.null(plottitle)) plottitle = paste(chromosome,": ",Start,"-",End,sep="")
	title(plottitle,cex.main=2)
	legend(legendloc,names(groupsplit),lty=c(1,1),col=c("blue","red"))
	if(skiplines){
		for(i in which(diff(pos)[xaxinds]!=1)[-length(which(diff(pos)[xaxinds]!=1))]){
		   abline(v=xaxinds[i],lwd=0.5,col="black")}
	}
	

	## panel B: plots of t-statistics
	plot(xaxinds,tstats[xaxinds],type="l",col="gray39",lwd=3,xaxt="n",xlim=c(plotlb,plotub),xlab="",ylab="")
	title(ylab="t statistic",mgp=c(1.5,1,0.2))
	abline(h=0,lty=2)
	if(skiplines){
		for(i in which(diff(pos)[xaxinds]!=1)[-length(which(diff(pos)[xaxinds]!=1))]){
		   abline(v=xaxinds[i],lwd=0.5,col="black")}
	}
	
	
	## panel C: annotated exons and state calls.
	# may want to fix this later to better handle overlapping annotated exons.
	#setup
	plot(xaxinds,rep(0,length(xaxinds)),type="n",ylim=c(-.5,.5),yaxt="n",ylab="",xlim=c(plotlb,plotub),xaxt="n")
	#exons
	exonInds = which(annotation$end >= lowerbound & annotation$start <= upperbound & annotation$chr == chromosome )
	if(length(exonInds)==0){warning("FYI: no annotated exons in specified region")}
	theExons = annotation[exonInds,]
	for(i in exonInds){
		exonend <- which(pos==annotation$end[i])
		exonstart <- which(pos==annotation$start[i])
		if(length(exonend)==0) {
			exonend <- which(pos<=annotation$end[i])[length(which(pos<=annotation$end[i]))]
			warning(paste("exon end position not in pos. Actual end: ",annotation$end[i],", marked end: ", pos[exonend],sep=""))
		}
		if(length(exonstart)==0){
			exonstart <- which(pos>=annotation$start[i])[1]
			warning(paste("exon start position not in pos. Actual start: ",annotation$start[i],", marked start: ", pos[exonstart],sep=""))
		}
		polygon(x=c(rep(exonstart,2),rep(exonend,2)), y=c(-.1,-.4,-.4,-.1), col="purple3")
	} # end plotting of exons
	#states
	for(k in 1:length(xaxinds)){
		if(getRegionObject$states.norle$state[xaxinds][k]==1) thecolor <- "gray"
		if(getRegionObject$states.norle$state[xaxinds][k]==2) thecolor <- "black"
		if(getRegionObject$states.norle$state[xaxinds][k]==3) thecolor <- "red"
		if(getRegionObject$states.norle$state[xaxinds][k]==4) thecolor <- "green"
		lines(rep(xaxinds[k],2),c(.1,.4),col=thecolor,lwd=2)
	} # end plotting of states
	axis(2,at=c(-.25,.25),labels=c("exons","states"),cex.axis=1.5)
	axsize = ifelse(is.null(plotfile),1.5,2)
	axis(1,at=xaxinds[seq(1,length(xaxinds),by=100)],labels=pos[xaxinds[seq(1,length(xaxinds),by=100)]],cex.axis=axsize)
	title(xlab="genomic position",outer=TRUE)
	if(skiplines){
		for(i in which(diff(pos)[xaxinds]!=1)[-length(which(diff(pos)[xaxinds]!=1))]){
		   abline(v=xaxinds[i],lwd=0.5,col="black")}
	}
	
	if(!is.null(plotfile)) dev.off()
}

	
