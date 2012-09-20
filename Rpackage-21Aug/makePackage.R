## create alpha version of R package
## AF 8/21/12

source("~/Documents/hopkins/research/_readlets/tornado/getAnnotation.R")
source("~/Documents/hopkins/research/_readlets/tornado/getExons.R")
source("~/Documents/hopkins/research/_readlets/tornado/getLimmaInput-parallel.R")
source("~/Documents/hopkins/research/_readlets/tornado/getParams.R")
source("~/Documents/hopkins/research/_readlets/tornado/getRegions.R")
source("~/Documents/hopkins/research/_readlets/tornado/getTstats.R")
source("~/Documents/hopkins/research/_readlets/tornado/locfdrFit.R")
source("~/Documents/hopkins/research/_readlets/tornado/makeDb.R")
source("~/Documents/hopkins/research/_readlets/tornado/paramHelpers.R")
source("~/Documents/hopkins/research/_readlets/tornado/plotfunctions.R")
source("~/Documents/hopkins/research/_readlets/tornado/sqldf-hacked.R")

package.skeleton(name="tornado")


getAnnotation = function (genome, tablename, genes = TRUE, verbose = TRUE) 
{
    require(GenomicFeatures)
    require(rtracklayer)
    a <- try(makeTranscriptDbFromUCSC(genome = genome, tablename = tablename))
    if (class(a) == "try-error") 
        stop("Problem accessing requested UCSC annotation - likely there is a problem with genome or tablename arguments. Use ucscGenomes() to see acceptable genomes; use supportedTables(genome) to see acceptable tablenames for your genome.")
    if (verbose) 
        show(a)
    if (genes) {
        if (verbose) 
            print("Labeling exons by gene...")
        grl <- exonsBy(a, "gene")
        if (length(grl) == 0) 
            stop("cannot list exons by gene using this table. please try another table or specify genes = FALSE to list exons by transcript.")
        datf <- IRanges:::as.data.frame(grl)
        names(datf)[1] <- "gene"
    }
    if (!genes) {
        if (verbose) 
            print("Labeling exons by transcript...")
        grl <- exonsBy(a, "tx")
        datf <- IRanges:::as.data.frame(grl)
        names(datf)[1] <- "transcript"
        nametable <- id2name(a, feature.type = "tx")
    }
    names(datf)[2] <- "chr"
    datf <- datf[, 1:7]
    if (genes) 
        return(datf)
    return(list(exon.table = datf, transcript.names = nametable))
}
