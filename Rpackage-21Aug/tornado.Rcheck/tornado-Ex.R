pkgname <- "tornado"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('tornado')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("getAnnotation")
### * getAnnotation

flush(stderr()); flush(stdout())

### Name: getAnnotation
### Title: download exon information for a given genome
### Aliases: getAnnotation

### ** Examples

mouse.exons <- getAnnotation("mm9","refGene")
head(mouse.exons)



cleanEx()
nameEx("getExons")
### * getExons

flush(stderr()); flush(stdout())

### Name: getExons
### Title: find closest exon(s) to a genomic region
### Aliases: getExons

### ** Examples

## not run:
exons = getAnnotation("hg19","knownGene")
theRegion = c("chr22", 18216902, 18218350)
getExons(theRegion, exons)
foo = getExons(theRegion, exons)
foo
foo$closestExons



cleanEx()
nameEx("getLimmaInput")
### * getLimmaInput

flush(stderr()); flush(stdout())

### Name: getLimmaInput
### Title: fit a linear model to each nucleotide
### Aliases: getLimmaInput

### ** Examples

## add example here when we have a vignette



cleanEx()
nameEx("last")
### * last

flush(stderr()); flush(stdout())

### Name: last
### Title: get the last element
### Aliases: last

### ** Examples

x = c(1:20)
last(x) #returns 20



cleanEx()
nameEx("supportedTables")
### * supportedTables

flush(stderr()); flush(stdout())

### Name: supportedTables
### Title: print list of supported (downloadable) tables for a given genome
### Aliases: supportedTables

### ** Examples

supportedTables("mm10")
mouse.exons = getAnnotation("mm10","refGene") #refGene appears in printed output of supportedTables("mm10").



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
