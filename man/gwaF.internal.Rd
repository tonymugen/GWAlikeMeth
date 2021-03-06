% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gwaF.internal}
\alias{gwaF.internal}
\title{GWA with fixed effects}
\usage{
gwaF.internal(yVec, kVec, xVec, snps, d, Ngen, nThr)
}
\arguments{
\item{yVec}{vectorized matrix of phenotypes}

\item{kVec}{vectorized relationship matrix}

\item{xVec}{vectorized fixed effect matrix}

\item{snps}{SNP matrix, SNPs as columns}

\item{d}{number of traits}

\item{Ngen}{number of genotypes}

\item{nThr}{number of threads}
}
\description{
Fits a random-effects model (with fixed effect covariates but no replicated measurement of genotypes) and does GWA on the provided SNPs. Operates on any number of traits at once, but treats them as independent. If the number of threads is set to zero, all available cores are used.
}
