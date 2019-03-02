/*
 * Copyright (c) 2018 Anthony J. Greenberg
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */


/// GWA on replicated data with a mixed model
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2017 Anthony J. Greenberg
 * \version 0.1
 *
 * R interface functions that perform GWA on replicated data. Fixed-effect covariates and missing SNP data are allowed. SNPs have to be coded as (0,1,2) with missing data marked as -9.
 * The implementation depends on C++-11. SNP regression is multi-threaded. Vectorized matrices can be passed directly from R, no trasition to row-major storage needed.
 * Multiple traits from the same experiment are treated at once, but the statistics are calculated independently (covariances effectively set to zero).
 *
 */

#include <vector>
#include <cmath>
#include <algorithm>
#include <thread>

#include <Rcpp.h>

#include "locMatrix.hpp"
#include "likeMeth.hpp"

//' Random effects fit
//'
//' Fits a random-effects model (with no fixed effect covariates other than the intercept). Operates on any number of traits at once, but treats them as independent.
//'
//' @param yVec vectorized matrix of phenotypes
//' @param kVec vectorized relationship matrix
//' @param repFac factor relating genotypes to replicates
//' @param d number of traits
//' @param Ngen number of genotypes
//' @export
//[[Rcpp::export(name="reFit")]]
Rcpp::List reFit(const std::vector<double> &yVec, const std::vector<double> &kVec, const std::vector<int32_t> &repFac, const int32_t &d, const int32_t &Ngen){
	if(d <= 0){
		Rcpp::stop("The number of traits must be positive");
	} else if (Ngen <= 0) {
		Rcpp::stop("The number of genotypes must be positive");
	}

	std::vector<size_t> fixedFac;
	for (auto &i : repFac) {
		if (i <= 0) {
			Rcpp::stop("Factor elements must be strictly positive");
		}
		fixedFac.push_back(static_cast<size_t>(i-1));
	}
	try {
		BayesicSpace::MixedModel model(yVec, kVec, fixedFac, d, Ngen, yVec.size()/d);

		BayesicSpace::Matrix U;    // matrix of random effects
		BayesicSpace::Matrix B;    // matrix of intercepts (no other fixed effects)
		std::vector<double> hSq;   // marker heritability
		std::vector<double> uVec;  // vectorized matrix of random effects
		std::vector<double> muVec; // vector of means
		model.ranef(U);
		U.vectorize(uVec);
		model.fixef(B);
		B.vectorize(muVec);
		model.hSq(hSq);
		return Rcpp::List::create(Rcpp::Named("ranef", uVec), Rcpp::Named("fixef", muVec), Rcpp::Named("hSq", hSq));
	} catch(std::string problem) {
		Rcpp::stop(problem);
	}

	return Rcpp::List::create(Rcpp::Named("error", "NaN"));
}

//' GWA with no fixed effect covariates
//'
//' Fits a random-effects model (with no fixed effect covariates other than the intercept) and does GWA on the provided SNPs. Operates on any number of traits at once, but treats them as independent.
//'
//' @param yVec vectorized matrix of phenotypes
//' @param kVec vectorized relationship matrix
//' @param repFac factor relating genotypes to replicates
//' @param snps SNP matrix, SNPs as columns
//' @param d number of traits
//' @param Ngen number of genotypes
//' @export
//[[Rcpp::export(name="gwa.internal")]]
Rcpp::List gwa(const std::vector<double> &yVec, const std::vector<double> &kVec, const std::vector<int32_t> &repFac, const std::vector<int32_t> &snps, const int32_t &d, const int32_t &Ngen){
	if(d <= 0){
		Rcpp::stop("The number of traits must be positive");
	} else if (Ngen <= 0) {
		Rcpp::stop("The number of genotypes must be positive");
	}

	std::vector<size_t> fixedFac;
	for (auto &i : repFac) {
		if (i <= 0) {
			Rcpp::stop("Factor elements must be strictly positive");
		}
		fixedFac.push_back(static_cast<size_t>(i-1));
	}
	try {
		std::vector<double> lPval;
		BayesicSpace::MixedModel model(yVec, kVec, fixedFac, d, Ngen, yVec.size()/d, &snps, -9, &lPval);

		BayesicSpace::Matrix U;    // matrix of random effects
		BayesicSpace::Matrix B;    // matrix of intercepts (no other fixed effects)
		std::vector<double> hSq;   // marker heritability
		std::vector<double> uVec;  // vectorized matrix of random effects
		std::vector<double> muVec; // vector of means
		model.ranef(U);
		U.vectorize(uVec);
		model.fixef(B);
		B.vectorize(muVec);
		model.hSq(hSq);
		model.gwa();
		return Rcpp::List::create(Rcpp::Named("ranef", uVec), Rcpp::Named("fixef", muVec), Rcpp::Named("hSq", hSq), Rcpp::Named("lPval", lPval));
	} catch(std::string problem) {
		Rcpp::stop(problem);
	}

	return Rcpp::List::create(Rcpp::Named("error", "NaN"));
}

//' GWA with FDR and no fixed effect covariates
//'
//' Fits a random-effects model (with no fixed effect covariates other than the intercept) and does GWA on the provided SNPs. Operates on any number of traits at once, but treats them as independent. Permutes the rows of the trait matrix to generate a null distribution of \f$ -\log_{10}p \f$ values. Uses this distribution to estimate per-SNP empirical false discovery rates.
//'
//' @param yVec vectorized matrix of phenotypes
//' @param kVec vectorized relationship matrix
//' @param repFac factor relating genotypes to replicates
//' @param snps SNP matrix, SNPs as columns
//' @param d number of traits
//' @param Ngen number of genotypes
//' @param nPer number of permutations
//' @export
//[[Rcpp::export(name="gwaFDR.internal")]]
Rcpp::List gwaFDR(const std::vector<double> &yVec, const std::vector<double> &kVec, const std::vector<int32_t> &repFac, const std::vector<int32_t> &snps, const int32_t &d, const int32_t &Ngen, const int32_t &nPer){
	if(d <= 0){
		Rcpp::stop("The number of traits must be positive");
	} else if (Ngen <= 0) {
		Rcpp::stop("The number of genotypes must be positive");
	} else if (nPer <= 0) {
		Rcpp::stop("The number of permutations must be positive");
	}

	std::vector<size_t> fixedFac;
	for (auto &i : repFac) {
		if (i <= 0) {
			Rcpp::stop("Factor elements must be strictly positive");
		}
		fixedFac.push_back(static_cast<size_t>(i-1));
	}
	try {
		std::vector<double> lPval;
		std::vector<double> fdr;
		BayesicSpace::MixedModel model(yVec, kVec, fixedFac, d, Ngen, yVec.size()/d, &snps, -9, &lPval);

		BayesicSpace::Matrix U;    // matrix of random effects
		BayesicSpace::Matrix B;    // matrix of intercepts (no other fixed effects)
		std::vector<double> hSq;   // marker heritability
		std::vector<double> uVec;  // vectorized matrix of random effects
		std::vector<double> muVec; // vector of means
		model.ranef(U);
		U.vectorize(uVec);
		model.fixef(B);
		B.vectorize(muVec);
		model.hSq(hSq);
		model.gwa(nPer, fdr);

		return Rcpp::List::create(Rcpp::Named("ranef", uVec), Rcpp::Named("fixef", muVec), Rcpp::Named("hSq", hSq), Rcpp::Named("lPval", lPval), Rcpp::Named("qVal", fdr));
	} catch(std::string problem) {
		Rcpp::stop(problem);
	}

	return Rcpp::List::create(Rcpp::Named("error", "NaN"));
}



