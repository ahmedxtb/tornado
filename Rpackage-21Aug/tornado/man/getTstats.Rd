\name{getTstats}
\alias{getTstats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ calculate moderated t statistic for each nucleotide
}
\description{
Less-complex version of eBayes in the limma package -- only gives output for one covariate of interest. The simplifying was done to conserve memory.
}
\usage{
getTstats(fit, trend = FALSE)
}
\arguments{
  \item{fit}{
list with elements \code{$sigma}, \code{$df.residual}, \code{$Amean} (if \code{trend=TRUE}), \code{$coefficients}, and \code{$stdev.unscaled} -- this can be an \code{lmFit} object or an object produced with \code{getLimmaInput}, namely, the \code{$ebobject} component of a \code{getLimmaInput} object.
}
  \item{trend}{
logical, should an intensity-trend be allowed for the prior variance? (see limma's eBayes). Default is that the prior variance is constant.
}
}
\value{
list with elements
\item{tt }{vector containing moderated t-statistics for each nucleotide}
\item{logfchange }{vector containing the log base 2 fold change in coverage between the two groups, as estimated by a linear model}
}
\references{
http://www.bioconductor.org/packages/2.10/bioc/vignettes/limma/inst/doc/usersguide.pdf
}
\author{
eBayes authored by Gordon Smyth, modifications by Alyssa Frazee
}

\seealso{
\code{\link{getLimmaInput}}
}
