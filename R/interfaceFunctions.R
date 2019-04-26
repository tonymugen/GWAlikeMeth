#
# Copyright (c) 2019 Anthony J. Greenberg
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
# BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE.
#

#' Mixed model fit
#'
#' Fits a mixed model using a relationship matrix and possible fixed effects (in addition to the intercept) and replication. Calculates breeding values and narrow-sense heritability. If the relationship matrix is built from genotypes, these are genome-estimated breeding values (GEBV) and marker heritability. Multiple traits are processed together, but treated independently in the model. Missing phenotype, fixed effect, or replication data are not supported. The Y and X parameters will be converted to one-column matrices if their type is not \code{matrix}. An attempt will be made to do convert K to a matrix with reasonable dimensions of a non-\code{matrix} type is passed. Do not include the intercept in the X matrix.
#'
#' @param Y matrix of phenotypes; phenotypes should be in columns.
#' @param K relationship matrix
#' @param X fixed effect matrix (optional, do not include the intercept)
#' @param repFactor factor relating genotypes to replicates (optional)
#' @return a list containing a matrix of random effects (named \code{ranef}), a matrix of fixed effects (including the intercept, named \code{fixef}), and a vector of heritability estimates (named \code{hSq})
#' @references
#' \insertRef{kang08}{GWAlikeMeth}
#' @export
mmFit <- function(Y, K, X=NULL, repFactor=NULL){
	if (!is.matrix(Y)) {
		Y <- as.matrix(Y)
	}
	d <- ncol(Y)
	# process the case with no replication
	if (is.null(repFactor)) {
		if (nrow(Y) != nrow(K)) {
			msg <- paste("Y and K must have the same number of rows if there is no replication; I have nrow(Y) = ", nrow(Y), " and nrow(K) = ", nrow(K), sep="")
			stop(msg)
		}
		if (!is.matrix(K)) {
			K <- matrix(K, nrow=nrow(Y))
			if (nrow(K) != ncol(K)) {
				msg <- paste("K must be square; I have nrow(K) = ", nrow(K), " and ncol(K) = ", ncol(K), sep="")
				stop(msg)
			}
		}
		if (is.null(X)) {
			res <- reFit(as.double(Y), as.double(K), d, nrow(K))
			res$ranef <- matrix(res$ranef, ncol=d)
			res$fixef <- matrix(res$fixef, ncol=d)
			return(res)
		} else {
			if (!is.matrix(X)) {
				X <- as.matrix(X)
				if (nrow(X) != nrow(Y)) {
					msg <- paste("X must have the same number of rows as Y when there is no replication; I have nrow(X) = ", nrow(X), " and nrow(Y) = ", nrow(Y), sep="")
					stop(msg)
				}
			}
			res <- reFitF(as.double(Y), as.double(K), as.double(X), d, nrow(K))
			res$ranef <- matrix(res$ranef, ncol=d)
			res$fixef <- matrix(res$fixef, ncol=d)
			return(res)
		}
	} else {
		# process the case with replication
		if (!is.factor(repFactor)) {
			repFactor <- as.factor(repFactor)
		}
		if ( nrow(Y) != length(repFactor) ) {
			msg <- paste("The length of the replication factor must be the same as the number of rows in Y; I have nrow(Y) = ", nrow(Y), " and length(repFactor) = ", length(repFactor), sep="")
			stop(msg)
		}
		if (!is.matrix(K)) {
			K <- matrix(K, nrow=nlevels(repFactor))
			if ( nrow(K) != ncol(K) ) {
				msg <- paste("Cannot coerce non-matrix K into a square matrix with correct dimensions; coerced K is ", nrow(K), "x", ncol(K), sep="")
				stop(msg)
			}
		}
		if ( nrow(K) != nlevels(repFactor) ) {
			msg <- paste("The number of rows in K must be the same as the number of levels in the replication factor; I have nrow(K) = ", nrow(K), " and nlevels(repFactor) = ", nlevels(repFactor), sep="")
			stop(msg)
		}
		#process the case with no fixed effect
		if (is.null(X)) {
			res <- reFitR(as.double(Y), as.double(K), as.integer(repFactor), d, nrow(K))
			res$ranef <- matrix(res$ranef, ncol=d)
			res$fixef <- matrix(res$fixef, ncol=d)
			return(res)
		} else {
			if (!is.matrix(X)) {
				X <- as.matrix(X)
			}
			if ( nrow(X) != nrow(Y) ) {
				if ( nrow(X) == nrow(K) ) { # expand X if it is on the genotype scale
					X <- X[repFactor,]
				} else {
					msg <- paste("The number of rows in X must be the same as either Y or K; I have nrow(X) = ", nrow(X), ", nrow(Y) = ", nrow(Y), ", and nrow(K) = ", nrow(K), sep="")
					stop(msg)
				}
			}
			res <- reFitRF(as.double(Y), as.double(K), as.integer(repFactor), as.double(X), d, nrow(K))
			res$ranef <- matrix(res$ranef, ncol=d)
			res$fixef <- matrix(res$fixef, ncol=d)
			return(res)
		}
	}
}

#' Single-marker genome-wide association
#'
#' Performs single-marker GWA with a mixed-model population structure correction. The algorithm is EMMAX. Optionally supports replicated genotypes and fixed-effect covariates. Empirical \eqn{q}-values can be calculated using permutations is desired. Since all loci are used in the permutation test, the number of samples need not be very large. I get good results with 15 to 25 permutations. Phenotype, relationship matrix, covariate, and replication specification is the same as in \code{\link{mmFit}}. The genotype matrix may have missing data (specified as \code{NA}). Genotypes should be coded as {0, 1, 2}, with 1 as the heterozygote. The lines should be in rows, with the number of rows equal to the number of lines.
#'
#' It is possible to specify the number of cores to be used. Multithreading is employed only for the SNP effect estimation. If you compiled R with multithreaded BLAS and/or LAPACK, it is safe to use it with this function. The linear algebra functions are run on single threads. The specified number of threads will be taken as a limit. If there are too few genotypes, not all of the suggested cores may be deployed. If the number of threads is \code{NULL}, 0, or not specified, the maximum number of cores will be used, provided there are enough SNPs to do so.
#'
#' @param Y matrix of phenotypes; phenotypes should be in columns.
#' @param K relationship matrix
#' @param X fixed effect matrix (optional, do not include the intercept)
#' @param repFactor factor relating genotypes to replicates (optional)
#' @param snps matrix of genotypes (\code{NA} for missing data, lines as rows, SNPs as columns)
#' @param nPerm number of permutations to calculate \eqn{q}-values
#' @param nThreads suggested number of threads; if \code{NULL} or 0, maximum number will be used
#' @return a list containing a matrix of random effects (named \code{ranef}), a matrix of fixed effects (including the intercept, named \code{fixef}), a vector of heritability estimates (named \code{hSq}), matrix of \eqn{-\log_{10}p} values (named \code{lPval}), and an optinal matrix of \eqn{q}-values (named \code{qVal}).
#' @references
#' \insertRef{kang08}{GWAlikeMeth}
#' \insertRef{kang10}{GWAlikeMeth}
#' @export
gwa <- function(Y, K, X=NULL, repFactor=NULL, snps, nPerm=NULL, nThreads=NULL){
	if(is.null(nThreads)){
		nThreads <- 0
	}
	snps[is.na(snps)] <- -9
	if(!is.matrix(snps)){
		snps <- as.matrix(snps)
	}
	if(!is.matrix(Y)){
		Y <- as.matrix(Y)
	}
	d <- ncol(Y)
	# the big bifurcation: test if permutations are needed
	if(is.null(nPerm) || (nPerm == 0)){
		if (is.null(repFactor)) {
			if (nrow(Y) != nrow(K)) {
				msg <- paste("Y and K must have the same number of rows if there is no replication; I have nrow(Y) = ", nrow(Y), " and nrow(K) = ", nrow(K), sep="")
				stop(msg)
			}
			if (!is.matrix(K)) {
				K <- matrix(K, nrow=nrow(Y))
				if (nrow(K) != ncol(K)) {
					msg <- paste("K must be square; I have nrow(K) = ", nrow(K), " and ncol(K) = ", ncol(K), sep="")
					stop(msg)
				}
				if (nrow(K) != nrow(snps)) {
					msg <- paste("The SNP matrix must have the same number of rows as K; I have nrow(K) = ", nrow(K), " and nrow(snps) = ", nrow(snps), sep="")
					stop(msg)
				}
			}
			if (is.null(X)) {
				res <- gwa.internal(as.double(Y), as.double(K), as.integer(snps), d, nrow(K), as.integer(nThreads))
				res$ranef <- matrix(res$ranef, ncol=d)
				res$fixef <- matrix(res$fixef, ncol=d)
				res$lPval <- matrix(res$lPval, ncol=d)
				return(res)
			} else {
				if (!is.matrix(X)) {
					X <- as.matrix(X)
					if (nrow(X) != nrow(Y)) {
						msg <- paste("X must have the same number of rows as Y when there is no replication; I have nrow(X) = ", nrow(X), " and nrow(Y) = ", nrow(Y), sep="")
						stop(msg)
					}
				}
				res <- gwaF.internal(as.double(Y), as.double(K), as.double(X), as.integer(snps), d, nrow(K), as.integer(nThreads))
				res$ranef <- matrix(res$ranef, ncol=d)
				res$fixef <- matrix(res$fixef, ncol=d)
				res$lPval <- matrix(res$lPval, ncol=d)
				return(res)
			}
		} else {
			# process the case with replication
			if (!is.factor(repFactor)) {
				repFactor <- as.factor(repFactor)
			}
			if ( nrow(Y) != length(repFactor) ) {
				msg <- paste("The length of the replication factor must be the same as the number of rows in Y; I have nrow(Y) = ", nrow(Y), " and length(repFactor) = ", length(repFactor), sep="")
				stop(msg)
			}
			if (!is.matrix(K)) {
				K <- matrix(K, nrow=nlevels(repFactor))
				if ( nrow(K) != ncol(K) ) {
					msg <- paste("Cannot coerce non-matrix K into a square matrix with correct dimensions; coerced K is ", nrow(K), "x", ncol(K), sep="")
					stop(msg)
				}
				if (nrow(K) != nrow(snps)) {
					msg <- paste("The SNP matrix must have the same number of rows as K; I have nrow(K) = ", nrow(K), " and nrow(snps) = ", nrow(snps), sep="")
					stop(msg)
				}
			}
			if ( nrow(K) != nlevels(repFactor) ) {
				msg <- paste("The number of rows in K must be the same as the number of levels in the replication factor; I have nrow(K) = ", nrow(K), " and nlevels(repFactor) = ", nlevels(repFactor), sep="")
				stop(msg)
			}
			#process the case with no fixed effect
			if (is.null(X)) {
				res <- gwaR.internal(as.double(Y), as.double(K), as.integer(repFactor), as.integer(snps), d, nrow(K), as.integer(nThreads))
				res$ranef <- matrix(res$ranef, ncol=d)
				res$fixef <- matrix(res$fixef, ncol=d)
				res$lPval <- matrix(res$lPval, ncol=d)
				return(res)
			} else {
				if (!is.matrix(X)) {
					X <- as.matrix(X)
				}
				if ( nrow(X) != nrow(Y) ) {
					if ( nrow(X) == nrow(K) ) { # expand X if it is on the genotype scale
						X <- X[repFactor,]
					} else {
						msg <- paste("The number of rows in X must be the same as either Y or K; I have nrow(X) = ", nrow(X), ", nrow(Y) = ", nrow(Y), ", and nrow(K) = ", nrow(K), sep="")
						stop(msg)
					}
				}
				res <- gwaRF.internal(as.double(Y), as.double(K), as.integer(repFactor), as.double(X), as.integer(snps), d, nrow(K), as.integer(nThreads))
				res$ranef <- matrix(res$ranef, ncol=d)
				res$fixef <- matrix(res$fixef, ncol=d)
				res$lPval <- matrix(res$lPval, ncol=d)
				return(res)
			}
		}
	} else { # doing the permutations
		if (is.null(repFactor)) {
			if (nrow(Y) != nrow(K)) {
				msg <- paste("Y and K must have the same number of rows if there is no replication; I have nrow(Y) = ", nrow(Y), " and nrow(K) = ", nrow(K), sep="")
				stop(msg)
			}
			if (!is.matrix(K)) {
				K <- matrix(K, nrow=nrow(Y))
				if (nrow(K) != ncol(K)) {
					msg <- paste("K must be square; I have nrow(K) = ", nrow(K), " and ncol(K) = ", ncol(K), sep="")
					stop(msg)
				}
				if (nrow(K) != nrow(snps)) {
					msg <- paste("The SNP matrix must have the same number of rows as K; I have nrow(K) = ", nrow(K), " and nrow(snps) = ", nrow(snps), sep="")
					stop(msg)
				}
			}
			if (is.null(X)) {
				res <- gwaFDR.internal(as.double(Y), as.double(K), as.integer(snps), d, nrow(K), as.integer(nPerm), as.integer(nThreads))
				res$ranef <- matrix(res$ranef, ncol=d)
				res$fixef <- matrix(res$fixef, ncol=d)
				res$lPval <- matrix(res$lPval, ncol=d)
				res$qVal  <- matrix(res$qVal, ncol=d)
				return(res)
			} else {
				if (!is.matrix(X)) {
					X <- as.matrix(X)
					if (nrow(X) != nrow(Y)) {
						msg <- paste("X must have the same number of rows as Y when there is no replication; I have nrow(X) = ", nrow(X), " and nrow(Y) = ", nrow(Y), sep="")
						stop(msg)
					}
				}
				res <- gwaFDRF.internal(as.double(Y), as.double(K), as.double(X), as.integer(snps), d, nrow(K), as.integer(nPerm), as.integer(nThreads))
				res$ranef <- matrix(res$ranef, ncol=d)
				res$fixef <- matrix(res$fixef, ncol=d)
				res$lPval <- matrix(res$lPval, ncol=d)
				res$qVal  <- matrix(res$qVal, ncol=d)
				return(res)
			}
		} else {
			# process the case with replication
			if (!is.factor(repFactor)) {
				repFactor <- as.factor(repFactor)
			}
			if ( nrow(Y) != length(repFactor) ) {
				msg <- paste("The length of the replication factor must be the same as the number of rows in Y; I have nrow(Y) = ", nrow(Y), " and length(repFactor) = ", length(repFactor), sep="")
				stop(msg)
			}
			if (!is.matrix(K)) {
				K <- matrix(K, nrow=nlevels(repFactor))
				if ( nrow(K) != ncol(K) ) {
					msg <- paste("Cannot coerce non-matrix K into a square matrix with correct dimensions; coerced K is ", nrow(K), "x", ncol(K), sep="")
					stop(msg)
				}
				if (nrow(K) != nrow(snps)) {
					msg <- paste("The SNP matrix must have the same number of rows as K; I have nrow(K) = ", nrow(K), " and nrow(snps) = ", nrow(snps), sep="")
					stop(msg)
				}
			}
			if ( nrow(K) != nlevels(repFactor) ) {
				msg <- paste("The number of rows in K must be the same as the number of levels in the replication factor; I have nrow(K) = ", nrow(K), " and nlevels(repFactor) = ", nlevels(repFactor), sep="")
				stop(msg)
			}
			#process the case with no fixed effect
			if (is.null(X)) {
				res <- gwaFDRR.internal(as.double(Y), as.double(K), as.integer(repFactor), as.integer(snps), d, nrow(K), as.integer(nPerm), as.integer(nThreads))
				res$ranef <- matrix(res$ranef, ncol=d)
				res$fixef <- matrix(res$fixef, ncol=d)
				res$lPval <- matrix(res$lPval, ncol=d)
				res$qVal  <- matrix(res$qVal, ncol=d)
				return(res)
			} else {
				if (!is.matrix(X)) {
					X <- as.matrix(X)
				}
				if ( nrow(X) != nrow(Y) ) {
					if ( nrow(X) == nrow(K) ) { # expand X if it is on the genotype scale
						X <- X[repFactor,]
					} else {
						msg <- paste("The number of rows in X must be the same as either Y or K; I have nrow(X) = ", nrow(X), ", nrow(Y) = ", nrow(Y), ", and nrow(K) = ", nrow(K), sep="")
						stop(msg)
					}
				}
				res <- gwaFDRRF.internal(as.double(Y), as.double(K), as.integer(repFactor), as.double(X), as.integer(snps), d, nrow(K), as.integer(nPerm), as.integer(nThreads))
				res$ranef <- matrix(res$ranef, ncol=d)
				res$fixef <- matrix(res$fixef, ncol=d)
				res$lPval <- matrix(res$lPval, ncol=d)
				res$qVal  <- matrix(res$qVal, ncol=d)
				return(res)
			}
		}

	}
}


