#' Empirical Bayes Statistics For Differential Expression
#'
#' Computes empirical Bayes moderated-t p-values relative to a minimum meaningful fold-change threshold.
#'
#' This is an internal function adapted from the \code{treat} function in the \code{limma} package. It is used to test for a log2 fold change value (Dam-fusion/Dam) greater than the specified threshold.
#'

daimTreat <- function(fit, lfc = log2(1.2), trend=FALSE, robust=FALSE, winsor.tail.p=c(0.05,0.1))
    #	Moderated t-statistics with threshold
    #	Davis McCarthy, Gordon Smyth
    #	25 July 2008.  Last revised 27 February 2014.
{
    #	Check fit
    if(!is(fit,"MArrayLM")) stop("fit must be an MArrayLM object")
    if(is.null(fit$coefficients)) stop("coefficients not found in fit object")
    if(is.null(fit$stdev.unscaled)) stop("stdev.unscaled not found in fit object")
    fit$lods <- NULL

    coefficients <- as.matrix(fit$coefficients)
    stdev.unscaled <- as.matrix(fit$stdev.unscaled)
    sigma <- fit$sigma
    df.residual <- fit$df.residual
    if (is.null(coefficients) || is.null(stdev.unscaled) || is.null(sigma) ||
        is.null(df.residual))
        stop("No data, or argument is not a valid lmFit object")
    if (all(df.residual == 0))
        stop("No residual degrees of freedom in linear model fits")
    if (all(!is.finite(sigma)))
        stop("No finite residual standard deviations")
    if(trend) {
        covariate <- fit$Amean
        if(is.null(covariate)) stop("Need Amean component in fit to estimate trend")
    } else {
        covariate <- NULL
    }
    sv <- limma::squeezeVar(sigma^2, df.residual, covariate=covariate, robust=robust, winsor.tail.p=winsor.tail.p)
    fit$df.prior <- sv$df.prior
    fit$s2.prior <- sv$var.prior
    fit$s2.post <- sv$var.post
    df.total <- df.residual + sv$df.prior
    df.pooled <- sum(df.residual,na.rm=TRUE)
    df.total <- pmin(df.total,df.pooled)
    fit$df.total <- df.total
    se <- stdev.unscaled*sqrt(fit$s2.post)
    tstat <- (coefficients-lfc)/se
    fit$t <- array(0,dim(coefficients),dimnames=dimnames(coefficients))
    fit$p.value <- pt(tstat, df=df.total, lower.tail=FALSE)
    tstat <- pmax(tstat,0)
    fit$t <- tstat
    fit$treat.lfc <- lfc
    fit
}
