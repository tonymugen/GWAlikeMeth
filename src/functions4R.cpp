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
//' Fits a random-effects model (with no fixed effect covariates other than the intercept and no replication). Operates on any number of traits at once, but treats them as independent.
//'
//' @param yVec vectorized matrix of phenotypes
//' @param kVec vectorized relationship matrix
//' @param d number of traits
//' @param Ngen number of genotypes
//' @export
//[[Rcpp::export(name="reFit")]]
Rcpp::List reFit(const std::vector<double> &yVec, const std::vector<double> &kVec, const int32_t &d, const int32_t &Ngen){
	if(d <= 0){
		Rcpp::stop("The number of traits must be positive");
	} else if (Ngen <= 0) {
		Rcpp::stop("The number of genotypes must be positive");
	}

	try {
		BayesicSpace::MixedModel model(yVec, kVec, d, Ngen);

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

//' Random effects fit with replication
//'
//' Fits a random-effects model (with no fixed effect covariates other than the intercept). Operates on any number of traits at once, but treats them as independent.
//'
//' @param yVec vectorized matrix of phenotypes
//' @param kVec vectorized relationship matrix
//' @param repFac factor relating genotypes to replicates
//' @param d number of traits
//' @param Ngen number of genotypes
//' @export
//[[Rcpp::export(name="reFitR")]]
Rcpp::List reFitR(const std::vector<double> &yVec, const std::vector<double> &kVec, const std::vector<int32_t> &repFac, const int32_t &d, const int32_t &Ngen){
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

//' Random effects fit with fixed effects but no replication
//'
//' Fits a random-effects model (with fixed effect covariates). Operates on any number of traits at once, but treats them as independent.
//'
//' @param yVec vectorized matrix of phenotypes
//' @param kVec vectorized relationship matrix
//' @param xVec vectorized matrix of fixed effects
//' @param d number of traits
//' @param Ngen number of genotypes
//' @export
//[[Rcpp::export(name="reFitF")]]
Rcpp::List reFitF(const std::vector<double> &yVec, const std::vector<double> &kVec, const std::vector<double> &xvec, const int32_t &d, const int32_t &Ngen){
	if(d <= 0){
		Rcpp::stop("The number of traits must be positive");
	} else if (Ngen <= 0) {
		Rcpp::stop("The number of genotypes must be positive");
	}

	try {
		BayesicSpace::MixedModel model(yVec, kVec, xvec, d, Ngen);

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

//' Random effects fit with fixed effects and replication
//'
//' Fits a random-effects model (with fixed effect covariates and replicated genotype measurements). Operates on any number of traits at once, but treats them as independent.
//'
//' @param yVec vectorized matrix of phenotypes
//' @param kVec vectorized relationship matrix
//' @param repFac factor relating genotypes to replicates
//' @param xVec vectorized matrix of fixed effects
//' @param d number of traits
//' @param Ngen number of genotypes
//' @export
//[[Rcpp::export(name="reFitRF")]]
Rcpp::List reFitRF(const std::vector<double> &yVec, const std::vector<double> &kVec, const std::vector<int32_t> &repFac, const std::vector<double> &xvec, const int32_t &d, const int32_t &Ngen){
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
		BayesicSpace::MixedModel model(yVec, kVec, fixedFac, xvec, d, Ngen, yVec.size()/d);

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

//' Simple GWA
//'
//' Fits a random-effects model (with no fixed effect covariates other than the intercept and no replicated measurement of genotypes) and does GWA on the provided SNPs. Operates on any number of traits at once, but treats them as independent. If the number of threads is set to zero, all available cores are used.
//'
//' @param yVec vectorized matrix of phenotypes
//' @param kVec vectorized relationship matrix
//' @param snps SNP matrix, SNPs as columns
//' @param d number of traits
//' @param Ngen number of genotypes
//' @param nThr number of threads
//' @export
//[[Rcpp::export(name="gwa.internal")]]
Rcpp::List gwa(const std::vector<double> &yVec, const std::vector<double> &kVec, const std::vector<int32_t> &snps, const int32_t &d, const int32_t &Ngen, const int32_t &nThr){
	if(d <= 0){
		Rcpp::stop("The number of traits must be positive");
	} else if (Ngen <= 0) {
		Rcpp::stop("The number of genotypes must be positive");
	} else if (nThr < 0) {
		Rcpp::stop("The number of threads must be non-negative");
	}

	try {
		std::vector<double> lPval;
		BayesicSpace::MixedModel model(yVec, kVec, d, Ngen, &snps, -9, &lPval);

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
		model.gwa(nThr);
		return Rcpp::List::create(Rcpp::Named("ranef", uVec), Rcpp::Named("fixef", muVec), Rcpp::Named("hSq", hSq), Rcpp::Named("lPval", lPval));
	} catch(std::string problem) {
		Rcpp::stop(problem);
	}

	return Rcpp::List::create(Rcpp::Named("error", "NaN"));
}

//' GWA with fixed effects
//'
//' Fits a random-effects model (with fixed effect covariates but no replicated measurement of genotypes) and does GWA on the provided SNPs. Operates on any number of traits at once, but treats them as independent. If the number of threads is set to zero, all available cores are used.
//'
//' @param yVec vectorized matrix of phenotypes
//' @param kVec vectorized relationship matrix
//' @param xVec vectorized fixed effect matrix
//' @param snps SNP matrix, SNPs as columns
//' @param d number of traits
//' @param Ngen number of genotypes
//' @param nThr number of threads
//' @export
//[[Rcpp::export(name="gwaF.internal")]]
Rcpp::List gwaF(const std::vector<double> &yVec, const std::vector<double> &kVec, const std::vector<double> &xVec, const std::vector<int32_t> &snps, const int32_t &d, const int32_t &Ngen, const int32_t &nThr){
	if(d <= 0){
		Rcpp::stop("The number of traits must be positive");
	} else if (Ngen <= 0) {
		Rcpp::stop("The number of genotypes must be positive");
	} else if (nThr < 0) {
		Rcpp::stop("The number of threads must be non-negative");
	}

	try {
		std::vector<double> lPval;
		BayesicSpace::MixedModel model(yVec, kVec, xVec, d, Ngen, &snps, -9, &lPval);

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
		model.gwa(nThr);
		return Rcpp::List::create(Rcpp::Named("ranef", uVec), Rcpp::Named("fixef", muVec), Rcpp::Named("hSq", hSq), Rcpp::Named("lPval", lPval));
	} catch(std::string problem) {
		Rcpp::stop(problem);
	}

	return Rcpp::List::create(Rcpp::Named("error", "NaN"));
}

//' GWA with replication
//'
//' Fits a random-effects model (with no fixed effect covariates other than the intercept but with genotypes measured in replicates) and does GWA on the provided SNPs. Operates on any number of traits at once, but treats them as independent. All available cores will be used if the number of threads is set to zero.
//'
//' @param yVec vectorized matrix of phenotypes
//' @param kVec vectorized relationship matrix
//' @param repFac factor relating genotypes to replicates
//' @param snps SNP matrix, SNPs as columns
//' @param d number of traits
//' @param Ngen number of genotypes
//' @param nThr number of threads
//' @export
//[[Rcpp::export(name="gwaR.internal")]]
Rcpp::List gwaR(const std::vector<double> &yVec, const std::vector<double> &kVec, const std::vector<int32_t> &repFac, const std::vector<int32_t> &snps, const int32_t &d, const int32_t &Ngen, const int32_t &nThr){
	if(d <= 0){
		Rcpp::stop("The number of traits must be positive");
	} else if (Ngen <= 0) {
		Rcpp::stop("The number of genotypes must be positive");
	} else if (nThr < 0){
		Rcpp::stop("The number of threads must be non-negative");
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
		model.gwa(nThr);
		return Rcpp::List::create(Rcpp::Named("ranef", uVec), Rcpp::Named("fixef", muVec), Rcpp::Named("hSq", hSq), Rcpp::Named("lPval", lPval));
	} catch(std::string problem) {
		Rcpp::stop(problem);
	}

	return Rcpp::List::create(Rcpp::Named("error", "NaN"));
}

//' GWA with replication and fixed effects
//'
//' Fits a random-effects model (with fixed effect covariates and genotypes measured in replicates) and does GWA on the provided SNPs. Operates on any number of traits at once, but treats them as independent. All available cores will be used if the number of threads is set to zero.
//'
//' @param yVec vectorized matrix of phenotypes
//' @param kVec vectorized relationship matrix
//' @param repFac factor relating genotypes to replicates
//' @param xVec vectorized relationship matrix
//' @param snps SNP matrix, SNPs as columns
//' @param d number of traits
//' @param Ngen number of genotypes
//' @param nThr number of threads
//' @export
//[[Rcpp::export(name="gwaRF.internal")]]
Rcpp::List gwaRF(const std::vector<double> &yVec, const std::vector<double> &kVec, const std::vector<int32_t> &repFac, const std::vector<double> xVec, const std::vector<int32_t> &snps, const int32_t &d, const int32_t &Ngen, const int32_t &nThr){
	if(d <= 0){
		Rcpp::stop("The number of traits must be positive");
	} else if (Ngen <= 0) {
		Rcpp::stop("The number of genotypes must be positive");
	} else if (nThr < 0){
		Rcpp::stop("The number of threads must be non-negative");
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
		BayesicSpace::MixedModel model(yVec, kVec, fixedFac, xVec, d, Ngen, yVec.size()/d, &snps, -9, &lPval);

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
		model.gwa(nThr);
		return Rcpp::List::create(Rcpp::Named("ranef", uVec), Rcpp::Named("fixef", muVec), Rcpp::Named("hSq", hSq), Rcpp::Named("lPval", lPval));
	} catch(std::string problem) {
		Rcpp::stop(problem);
	}

	return Rcpp::List::create(Rcpp::Named("error", "NaN"));
}

//' Simple GWA with FDR
//'
//' Fits a random-effects model (with no fixed effect covariates other than the intercept and no replication) and does GWA on the provided SNPs. Operates on any number of traits at once, but treats them as independent. Permutes the rows of the trait matrix to generate a null distribution of \f$ -\log_{10}p \f$ values. Uses this distribution to estimate per-SNP empirical false discovery rates. If the number of threads is set to 0, the number is picked automatically.
//'
//' @param yVec vectorized matrix of phenotypes
//' @param kVec vectorized relationship matrix
//' @param snps SNP matrix, SNPs as columns
//' @param d number of traits
//' @param Ngen number of genotypes
//' @param nPer number of permutations
//' @param nThr number of threads
//' @export
//[[Rcpp::export(name="gwaFDR.internal")]]
Rcpp::List gwaFDR(const std::vector<double> &yVec, const std::vector<double> &kVec, const std::vector<int32_t> &snps, const int32_t &d, const int32_t &Ngen, const int32_t &nPer, const int32_t &nThr){
	if(d <= 0){
		Rcpp::stop("The number of traits must be positive");
	} else if (Ngen <= 0) {
		Rcpp::stop("The number of genotypes must be positive");
	} else if (nPer <= 0) {
		Rcpp::stop("The number of permutations must be positive");
	} else if (nThr < 0) {
		Rcpp::stop("The number of threads must be non-negative");
	}

	try {
		std::vector<double> lPval;
		std::vector<double> fdr;
		BayesicSpace::MixedModel model(yVec, kVec, d, Ngen, &snps, -9, &lPval);

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
		model.gwa(nPer, nThr, fdr);

		return Rcpp::List::create(Rcpp::Named("ranef", uVec), Rcpp::Named("fixef", muVec), Rcpp::Named("hSq", hSq), Rcpp::Named("lPval", lPval), Rcpp::Named("qVal", fdr));
	} catch(std::string problem) {
		Rcpp::stop(problem);
	}

	return Rcpp::List::create(Rcpp::Named("error", "NaN"));
}

//' GWA with FDR and replication
//'
//' Fits a random-effects model (with no fixed effect covariates other than the intercept) and does GWA on the provided SNPs. Operates on any number of traits at once, but treats them as independent. Permutes the rows of the trait matrix to generate a null distribution of \f$ -\log_{10}p \f$ values. Uses this distribution to estimate per-SNP empirical false discovery rates. The number of threads is set automatically if the number provided is 0.
//'
//' @param yVec vectorized matrix of phenotypes
//' @param kVec vectorized relationship matrix
//' @param repFac factor relating genotypes to replicates
//' @param snps SNP matrix, SNPs as columns
//' @param d number of traits
//' @param Ngen number of genotypes
//' @param nPer number of permutations
//' @param nThr number of threads
//' @export
//[[Rcpp::export(name="gwaFDRR.internal")]]
Rcpp::List gwaFDRR(const std::vector<double> &yVec, const std::vector<double> &kVec, const std::vector<int32_t> &repFac, const std::vector<int32_t> &snps, const int32_t &d, const int32_t &Ngen, const int32_t &nPer, const int32_t &nThr){
	if(d <= 0){
		Rcpp::stop("The number of traits must be positive");
	} else if (Ngen <= 0) {
		Rcpp::stop("The number of genotypes must be positive");
	} else if (nPer <= 0) {
		Rcpp::stop("The number of permutations must be positive");
	} else if (nThr < 0) {
		Rcpp::stop("The number of threads must be positive");
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
		model.gwa(nPer, nThr, fdr);

		return Rcpp::List::create(Rcpp::Named("ranef", uVec), Rcpp::Named("fixef", muVec), Rcpp::Named("hSq", hSq), Rcpp::Named("lPval", lPval), Rcpp::Named("qVal", fdr));
	} catch(std::string problem) {
		Rcpp::stop(problem);
	}

	return Rcpp::List::create(Rcpp::Named("error", "NaN"));
}

//' GWA with FDR and fixed effects
//'
//' Fits a random-effects model (with no fixed effect covariates other than the intercept) and does GWA on the provided SNPs. Operates on any number of traits at once, but treats them as independent. Permutes the rows of the trait matrix to generate a null distribution of \f$ -\log_{10}p \f$ values. Uses this distribution to estimate per-SNP empirical false discovery rates. The number of threads is set automatically if the number provided is 0.
//'
//' @param yVec vectorized matrix of phenotypes
//' @param kVec vectorized relationship matrix
//' @param xVec vectorized matrix of fixed effects
//' @param snps SNP matrix, SNPs as columns
//' @param d number of traits
//' @param Ngen number of genotypes
//' @param nPer number of permutations
//' @param nThr number of threads
//' @export
//[[Rcpp::export(name="gwaFDRF.internal")]]
Rcpp::List gwaFDRF(const std::vector<double> &yVec, const std::vector<double> &kVec, const std::vector<double> &xVec, const std::vector<int32_t> &snps, const int32_t &d, const int32_t &Ngen, const int32_t &nPer, const int32_t &nThr){
	if(d <= 0){
		Rcpp::stop("The number of traits must be positive");
	} else if (Ngen <= 0) {
		Rcpp::stop("The number of genotypes must be positive");
	} else if (nPer <= 0) {
		Rcpp::stop("The number of permutations must be positive");
	} else if (nThr < 0) {
		Rcpp::stop("The number of threads must be positive");
	}

	try {
		std::vector<double> lPval;
		std::vector<double> fdr;
		BayesicSpace::MixedModel model(yVec, kVec, xVec, d, Ngen, &snps, -9, &lPval);

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
		model.gwa(nPer, nThr, fdr);

		return Rcpp::List::create(Rcpp::Named("ranef", uVec), Rcpp::Named("fixef", muVec), Rcpp::Named("hSq", hSq), Rcpp::Named("lPval", lPval), Rcpp::Named("qVal", fdr));
	} catch(std::string problem) {
		Rcpp::stop(problem);
	}

	return Rcpp::List::create(Rcpp::Named("error", "NaN"));
}

//' GWA with FDR, replication, and fixed effects
//'
//' Fits a random-effects model (with no fixed effect covariates other than the intercept) and does GWA on the provided SNPs. Operates on any number of traits at once, but treats them as independent. Permutes the rows of the trait matrix to generate a null distribution of \f$ -\log_{10}p \f$ values. Uses this distribution to estimate per-SNP empirical false discovery rates. The number of threads is set automatically if the number provided is 0.
//'
//' @param yVec vectorized matrix of phenotypes
//' @param kVec vectorized relationship matrix
//' @param repFac factor relating genotypes to replicates
//' @param xVec vectorized matrix of fixed effects
//' @param snps SNP matrix, SNPs as columns
//' @param d number of traits
//' @param Ngen number of genotypes
//' @param nPer number of permutations
//' @param nThr number of threads
//' @export
//[[Rcpp::export(name="gwaFDRRF.internal")]]
Rcpp::List gwaFDRRF(const std::vector<double> &yVec, const std::vector<double> &kVec, const std::vector<int32_t> &repFac, const std::vector<double> &xVec, const std::vector<int32_t> &snps, const int32_t &d, const int32_t &Ngen, const int32_t &nPer, const int32_t &nThr){
	if(d <= 0){
		Rcpp::stop("The number of traits must be positive");
	} else if (Ngen <= 0) {
		Rcpp::stop("The number of genotypes must be positive");
	} else if (nPer <= 0) {
		Rcpp::stop("The number of permutations must be positive");
	} else if (nThr < 0) {
		Rcpp::stop("The number of threads must be positive");
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
		BayesicSpace::MixedModel model(yVec, kVec, fixedFac, xVec, d, Ngen, yVec.size()/d, &snps, -9, &lPval);

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
		model.gwa(nPer, nThr, fdr);

		return Rcpp::List::create(Rcpp::Named("ranef", uVec), Rcpp::Named("fixef", muVec), Rcpp::Named("hSq", hSq), Rcpp::Named("lPval", lPval), Rcpp::Named("qVal", fdr));
	} catch(std::string problem) {
		Rcpp::stop(problem);
	}

	return Rcpp::List::create(Rcpp::Named("error", "NaN"));
}



