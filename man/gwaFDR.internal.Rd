% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gwaFDR.internal}
\alias{gwaFDR.internal}
\title{Simple GWA with FDR}
\usage{
gwaFDR.internal(yVec, kVec, snps, d, Ngen, nPer, nThr)
}
\arguments{
\item{yVec}{vectorized matrix of phenotypes}

\item{kVec}{vectorized relationship matrix}

\item{snps}{SNP matrix, SNPs as columns}

\item{d}{number of traits}

\item{Ngen}{number of genotypes}

\item{nPer}{number of permutations}

\item{nThr}{number of threads}
}
\description{
Fits a random-effects model (with no fixed effect covariates other than the intercept and no replication) and does GWA on the provided SNPs. Operates on any number of traits at once, but treats them as independent. Permutes the rows of the trait matrix to generate a null distribution of \eqn{-\log_{10}p} values. Uses this distribution to estimate per-SNP empirical false discovery rates. If the number of threads is set to 0, the number is picked automatically.
}
