## 7/23/12 AF
## output of getLimmaInput --> moderated t-statistics
## idea: something like fit = getLimmaInput(); modt = getTstats(fit)$tt
## most of this function is directly copied out of eBayes() from limma - just kept the parts we needed

getTstats = function(fit, trend = FALSE){
	# get d0 and s02 (prior parameters)
	sg2 <- fit$sigma^2
	sg2 <- pmax(sg2,1e-05*median(sg2)) 
	zg <- log(sg2)
	eg <- zg - digamma(fit$df.residual/2)+log(fit$df.residual/2)
	rm(zg);gc();gc()
	if(!trend){
      ebar <- mean(eg)
      G <- length(fit$sigma)
      d0 <- 2*trigammaInverse(mean((eg-ebar)^2*(G/(G-1)) - trigamma(fit$df.residual/2)))
      s02 <- exp(ebar+digamma(d0/2)-log(d0/2))
  	}
	if(trend){
      require(splines)
      design <- try(ns(fit$Amean, df = 4, intercept = TRUE), silent = TRUE)
      if (is(design, "try-error")) stop("Problem with Amean covariate; perhaps too few distinct values")
      fit1 <- lm.fit(design, eg)
      ebar <- fit1$fitted
      evar <- mean(fit1$residuals[-(1:fit1$rank)]^2)
      rm(fit1);gc();gc()
      evar <- evar - mean(trigamma(fit$df.residual/2))
      d0 <- 2*trigammaInverse(evar)
      s02 <- exp(ebar+digamma(d0/2)-log(d0/2))
      rm(ebar,fit1,design);gc();gc()
  	}

	sgtilde <- sqrt((d0*s02+fit$df.residual*sg2)/(d0+fit$df.residual))
	tt <- fit$coefficients/(sgtilde*fit$stdev.unscaled)
	logfchange = fit$coefficients
	return(list(tt=tt,logfchange=logfchange))
}
