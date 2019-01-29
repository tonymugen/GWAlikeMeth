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
 * R interface functions that perform GWA on replcated data. Fixed-effect covariates and missing SNP data are allowed. SNPs have to be coded as (0,1,2) with missing data marked as -9.
 * The implementation depends on C++-11. SNP regression is multi-threaded. Vectorized matrices can be passed directly from R, no trasition to row-major storage needed. 
 * Multiple traits from the same experiment are treated at once, but the statistics are calculated independently (covariances effectively set to zero).
 *
 * To compile for linking from R, use
 *
 *	R CMD SHLIB repMMgwa.cpp locMatrix.cpp utilities.cpp likeMeth.cpp
 *
 * Make sure that C++-11 is enabled, for example by including a Makevars file with a line CXX_STD = CXX11 in the compilation directory.
 */

#include <vector>
#include <cmath>
#include <algorithm>
#include <thread>

#include "locMatrix.hpp"
#include "likeMeth.hpp"
#include <R.h> // has to be declared here because of namespace conflicts


using std::vector;
using std::fill;
using std::thread;
using locMatrix::Matrix;
typedef locMatrix::Matrix Matrix;


/** GWA with no fixed effects
 *
 * No fixed effects other than the intercept, which does not have to be provided. Matrices are passed in a vectorized form from R.
 *
 * \param[in] Yvec phenotypes (\f$N \times d \f$)
 * \param[in] Kvec kinship matrix (\f$N_{ln} \times N_{ln} \f$)
 * \param[in] Zvec design matrix (can be identity; \f$N \times N_{ln} \f$)
 * \param[in] SNPvec SNP table with SNPs as columns (\f$N_{ln} \times N_{snp} \f$)
 * \param[in] N number of data points
 * \param[in] Nln number of lines
 * \param[in] d number of traits
 * \param[in] Nsnp number of SNPs
 * \param[out] bHat estimated intercepts (\f$1 \times d \f$)
 * \param[out] uHat GEBV (\f$N_{ln} \times d \f$)
 * \param[out] lPval \f$-\log_{10}p\f$ (\f$N_{snp} \times d\f$)
 *
 *
 */
extern "C" {
	void gwa(const double *Yvec, const double *Kvec, const double *Zvec, const int *SNPvec, const int *N, const int *Nln, const int *d, const int *Nsnp, double *bHat, double *uHat, double *lPval){
		Matrix Y(Yvec, *N, *d);
		Matrix Z(Zvec, *N, *Nln);
		Matrix K(Kvec, *Nln, *Nln);
		
		Matrix u;
		vector<double> mu;
		
		solveMM(Y, K, Z, u, mu);
		
		// calculate within-line means of Y
		Y.premultZt(Z);
		vector<double> nLn; // number of replicates per line
		Z.colSums(nLn);
		Y.colDivide(nLn);
		
		Y = Y - u;
		Y.rowSub(mu);
		
		memcpy(bHat, mu.data(), (*d)*sizeof(double));
		mu.resize(0);
		nLn.resize(0);
		for (size_t jTrt = 0; jTrt < (*d); jTrt++) {
			for (size_t iLn = 0; iLn < (*Nln); iLn++) {
				uHat[(*Nln)*jTrt + iLn] = u.getElem(iLn, jTrt);
			}
		}
		u.resize(0, 0);
		Z.resize(0, 0);
		K.resize(0, 0);
		
		// determine the number of theards to use. There typically are two threads per core. It is worth doing no fewer than 200 SNPs per process.
		unsigned long nCores        = thread::hardware_concurrency()/2;
		unsigned long minBlockSize  = 200;
		unsigned long maxNumThreads = (*Nsnp)/minBlockSize;
		unsigned long nSNPsDone     = 0;

		if ((maxNumThreads <= 1) || (nCores == 1)) {
			for (size_t iSNP = 0; iSNP < (*Nsnp); iSNP++) {
				snpReg(Y, SNPvec + iSNP*(*Nln), -9, lPval + iSNP*(*d));
			}
		} else{
			unsigned long nThr = (maxNumThreads <= nCores ? maxNumThreads : nCores);
			unsigned long blockSize = (*Nsnp)/nThr;
			vector<thread> threads(nThr - 1);     // the main thread will do one block
			for (auto thrIt = threads.begin(); thrIt != threads.end(); ++thrIt) {
				(*thrIt) = thread{SNPblock(Y, SNPvec, nSNPsDone, blockSize, lPval)}; // has to be done on the fly because I deleted the copy constructor
				nSNPsDone += blockSize;
			}
			if (nSNPsDone < (*Nsnp)) { // if there are remaining SNPs to do, process the rest
				SNPblock block(Y, SNPvec, nSNPsDone, (*Nsnp) - nSNPsDone, lPval);
				block();
			}
			
			for (auto thrIt = threads.begin(); thrIt != threads.end(); ++thrIt) {
				if (thrIt->joinable()) {
					thrIt->join();
				}
			}
			
		}
		
		
	}
	
}

/** GWA with fixed effects
 *
 * Fixed effect covariates are included. They should not include the intercept since that is included internally. The number of rows in the covariate matrix can be either the same as the total number of data points (i.e. rows in \f$Y\f$) or the same as the number of lines (i.e. the number of rows in \f$K\f$ or columns in \f$Z\f$).
 * Matrices are passed in a vectorized form from R.
 *
 * \param[in] Yvec phenotypes (\f$N \times d \f$)
 * \param[in] Kvec kinship matrix (\f$N_{ln} \times N_{ln} \f$)
 * \param[in] Zvec design matrix (can be identity; \f$N \times N_{ln} \f$)
 * \param[in] Xvec covariate matrix  (\f$N \times N_{cvc}\f$ or \f$N_{ln} \times N_{cvc} \f$)
 * \param[in] SNPvec SNP table with SNPs as columns (\f$N_{ln} \times N_{snp} \f$)
 * \param[in] N number of data points
 * \param[in] Ncvr number of rows in the covariate matrix
 * \param[in] Ncvc number of columns in the covariate matrix
 * \param[in] Nln number of lines
 * \param[in] d number of traits
 * \param[in] Nsnp number of SNPs
 * \param[out] bHat estimated fixed effects (must have room for intercepts: \f$N_{cvc}+1 \times d \f$)
 * \param[out] uHat GEBV (\f$N_{ln} \times d \f$)
 * \param[out] lPval \f$-\log_{10}p\f$ (\f$N_{snp} \times d\f$)
 *
 *
 */
extern "C" {
	void gwaCv(const double *Yvec, const double *Kvec, const double *Zvec, const double *Xvec, const int *SNPvec, const int *N, const int *Ncvr, const int *Ncvc, const int *Nln, const int *d, const int *Nsnp, double *bHat, double *uHat, double *lPval){
		// check covariate dimensions first
		if ( ((*Ncvr) != (*Nln)) && ((*Ncvr) != (*N)) ) {
			fill(bHat, bHat + (*d)*(*Ncvc + 1), nan(""));
			fill(uHat, uHat + (*d)*(*Nln), nan(""));
			fill(lPval, lPval + (*d)*(*Nsnp), 0.0);
			return;
		}
		
		Matrix Y(Yvec, *N, *d);
		Matrix Z(Zvec, *N, *Nln);
		Matrix K(Kvec, *Nln, *Nln);
		Matrix X(Xvec, *Ncvr, *Ncvc);
		Matrix u;
		Matrix Beta;
		
		solveMM(Y, K, Z, X, u, Beta);
		
		// calculate within-line means of Y
		Y.premultZt(Z);
		vector<double> nLn; // number of replicates per line
		Z.colSums(nLn);
		Y.colDivide(nLn);
		
		Matrix XwIcpt(1.0, X.getNrows(), 1);
		XwIcpt.appendCol(X);
		Beta.gemm(false, 1.0, XwIcpt, false, 1.0, u); // u = u + Xbeta
		Y = Y - u;
		
		for (size_t jTrt = 0; jTrt < (*d); jTrt++) {
			for (size_t iLn = 0; iLn < (*Nln); iLn++) {
				uHat[(*Nln)*jTrt + iLn] = u.getElem(iLn, jTrt);
			}
			for (size_t iCv = 0; iCv < (*Ncvc) + 1; iCv++) {
				bHat[(*Ncvc+1)*jTrt + iCv] = Beta.getElem(iCv, jTrt);
			}
		}

		// free memory from variables I no longer need
		nLn.resize(0);
		XwIcpt.resize(0, 0);
		u.resize(0, 0);
		Beta.resize(0, 0);
		Z.resize(0, 0);
		K.resize(0, 0);
		X.resize(0, 0);
		
		// determine the number of theards to use. There typically are two threads per core. It is worth doing no fewer than 200 SNPs per process.
		unsigned long nCores        = thread::hardware_concurrency()/2;
		unsigned long minBlockSize  = 200;
		unsigned long maxNumThreads = (*Nsnp)/minBlockSize;
		unsigned long nSNPsDone     = 0;
		
		if ((maxNumThreads <= 1) || (nCores == 1)) {
			for (size_t iSNP = 0; iSNP < (*Nsnp); iSNP++) {
				snpReg(Y, SNPvec + iSNP*(*Nln), -9, lPval + iSNP*(*d));
			}
		} else{
			unsigned long nThr = (maxNumThreads <= nCores ? maxNumThreads : nCores);
			unsigned long blockSize = (*Nsnp)/nThr;
			vector<thread> threads(nThr - 1);     // the main thread will do one block
			for (auto thrIt = threads.begin(); thrIt != threads.end(); ++thrIt) {
				(*thrIt) = thread{SNPblock(Y, SNPvec, nSNPsDone, blockSize, lPval)}; // has to be done on the fly because I deleted the copy constructor
				nSNPsDone += blockSize;
			}
			if (nSNPsDone < (*Nsnp)) { // if there are remaining SNPs to do, process the rest
				SNPblock block(Y, SNPvec, nSNPsDone, (*Nsnp) - nSNPsDone, lPval);
				block();
			}
			
			for (auto thrIt = threads.begin(); thrIt != threads.end(); ++thrIt) {
				if (thrIt->joinable()) {
					thrIt->join();
				}
			}
			
		}
		
		
	}
	
}


