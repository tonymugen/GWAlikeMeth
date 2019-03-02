/*
 * Copyright (c) 2019 Anthony J. Greenberg
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


/// Likelihood methods for quantitative genetics
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2019 Anthony J. Greenberg
 * \version 0.1
 *
 * This is the file containing function implementations.
 *
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <thread>

#include "likeMeth.hpp"
#include "utilities.hpp"

using std::cerr;
using std::endl;
using std::vector;
using std::thread;

using namespace BayesicSpace;

double EmmREML::operator()(const double &delta){

	double frac = 0.0; // the log(Sum(eta^2/(lambda + delta))) component
	double llam = 0.0; // the Sum(log(lambda + delta))
	double ld;
	for (size_t s = 0; s < etaSq_->getNrows(); s++) {
		ld    = (*lambda_)[s] + delta;
		frac += etaSq_->getElem(s, jCol_)/ld;
		llam += log(ld);
	}
	// final likelihood
	double llik = -static_cast<double>(lambda_->size())*log(frac) - llam;
	return llik;

}

// MixedModel methods

MixedModel::MixedModel(const vector<double> &yvec, const vector<double> &kvec, const vector<size_t> &repFac, const size_t &d, const size_t &Ngen, const size_t &N) : Y_{Matrix(yvec, N, d)}, K_{Matrix(kvec, Ngen, Ngen)}, delta_{vector<double>(d, 0.0)} {

	if(repFac.size() != N){
		throw("ERROR: factor length not equal to the number of data points in MixedModel constructor");
	}
	Z_ = Matrix(0.0, N, Ngen);
	for (size_t i = 0; i < repFac.size(); ++i) {
		Z_.setElem(i, repFac[i], 1.0);
	}

	Matrix ZKZt(K_);
	double invN = -1.0/static_cast<double>(Z_.getNrows());
	Matrix S(invN, Z_.getNrows(), Z_.getNrows());
	for (size_t i = 0; i < Z_.getNrows(); i++) {
		S.setElem(i, i, 1.0 + invN);
	}

	// make ZKZ'; it is a symmetric matrix
	ZKZt.premultZ(Z_);
	// will need ZK for calculating BLUPs; save it
	Matrix ZK(ZKZt);
	ZKZt.postmultZt(Z_);
	double offset = log(static_cast<double>(ZKZt.getNrows()));
	for (size_t ii = 0; ii < ZKZt.getNrows(); ii++) {
		double tmp = ZKZt.getElem(ii, ii) + offset;
		ZKZt.setElem(ii, ii, tmp);
	}

	Matrix SHS; // SHS matrix (equation 5 of Kang et al.)
	S.symm('l', 'r', 1.0, ZKZt, 0.0, SHS);
	SHS.symm('l', 'r', 1.0, S, 0.0, u_); // u is now SHS

	vector<double> lam;
	S.resize(1, 1);
	u_.eigen('l', u_.getNcols() - 1, S, lam); // S now U(SHS)

	for (auto lamIt = lam.begin(); lamIt != lam.end(); ++lamIt) {
		if ((*lamIt) - offset <= 1e-10) {
			(*lamIt) = 0.0;
		} else {
			(*lamIt) -= offset;
		}

	}

	Matrix Eta;
	Y_.gemm(true, 1.0, S, false, 0.0, Eta); // U'Y

	for (size_t iRow = 0; iRow < Eta.getNrows(); iRow++) {
		for (size_t jCol = 0; jCol < Eta.getNcols(); jCol++) {
			Eta.setElem(iRow, jCol, pow2(Eta.getElem(iRow, jCol)));
		}
	}

	EmmREML reml(Eta, lam, 0);
	SHS = Y_; // replace SHS with Y. SHS will be modified in the loop
	vector<double> XtH(ZKZt.getNrows());
	for (size_t ii = 0; ii < ZKZt.getNrows(); ii++) {
		ZKZt.setElem(ii, ii, ZKZt.getElem(ii, ii) - offset);
	}

	beta_ = Matrix(1, d);
	for (size_t jCol = 0; jCol < Y_.getNcols(); jCol++) {
		reml.setColID(jCol);
		double fMax = 0.0;
		maximizer(reml, 1e-6, delta_[jCol], fMax);

		u_ = ZKZt; // replacing Zu with ZKZ'
		for (size_t ii = 0; ii < u_.getNrows(); ii++) {
			u_.setElem(ii, ii, ZKZt.getElem(ii, ii) + delta_[jCol]);
		}
		u_.chol();
		u_.cholInv();         // now H^{-1}
		u_.colSums(XtH);      // X'H^{-1}, since X is a column of 1s (intercept)
		double XtHX    = 0.0; // X'H^{-1}X
		double betaLoc = 0.0;
		for (size_t iRow = 0; iRow < Y_.getNrows(); iRow++) {
			XtHX    += XtH[iRow];
			betaLoc += Y_.getElem(iRow, jCol) * XtH[iRow];
		}

		betaLoc /= XtHX;
		beta_.setElem(0, jCol, betaLoc);

		SHS.colSub(betaLoc, jCol);
		u_.gemc(false, 1.0, SHS, jCol, 0.0, XtH);
		SHS.setCol(jCol, XtH);

	}
	SHS.gemm(true, 1.0, ZK, false, 0.0, u_); // replacing Zu with matrix of u
}

MixedModel::MixedModel(const vector<double> &yvec, const vector<double> &kvec, const vector<size_t> &repFac, const vector<double> &xvec, const size_t &d, const size_t &Ngen, const size_t &N) : Y_{Matrix(yvec, N, d)}, K_{Matrix(kvec, Ngen, Ngen)}, X_{Matrix(xvec, N, xvec.size()/N)}, delta_{vector<double>(d, 0.0)} {

	if (X_.getNcols() >= X_.getNrows()) {
		throw string("ERROR: covariate matrix X_ in MixedModel constructor is not full rank");
	}
	if ( (beta_.getNcols() != Y_.getNcols()) || (beta_.getNrows() != (X_.getNcols() + 1)) ) {
		beta_.resize(X_.getNcols() + 1, Y_.getNcols());
	}

	if(repFac.size() != N){
		throw("ERROR: factor length not equal to the number of data points in MixedModel constructor");
	}
	Z_ = Matrix(0.0, N, Ngen);
	for (size_t i = 0; i < repFac.size(); ++i) {
		Z_.setElem(i, repFac[i], 1.0);
	}
	// add the intercept column to X
	Matrix XwIcpt(1.0, X_.getNrows(), 1);
	XwIcpt.appendCol(X_);
	if (XwIcpt.getNrows() == K_.getNrows()) {
		XwIcpt.premultZ(Z_);
	} else if (XwIcpt.getNrows() != Y_.getNrows()) {
		throw string("ERROR: number of rows in X not equal to number of lines or number of data points");
	}
	Matrix S;
	XwIcpt.syrk('l', 1.0, 0.0, S);              // S = X'X
	S.chol();
	S.cholInv();                                // S = (X'X)^-1
	Matrix SHS;
	XwIcpt.symm('l', 'r', 1.0, S, 0.0, SHS);    // SHS = X(X'X)^-1
	XwIcpt.gemm(false, 1.0, SHS, true, 0.0, S); // S = SHS X' = X(X'X)^-1X'

	// I - S; lower triangle only
	for (size_t jCol = 0; jCol < S.getNcols(); jCol++) {
		S.setElem(jCol, jCol, 1.0 - S.getElem(jCol, jCol));
		for (size_t iRow = jCol+1; iRow < S.getNrows(); iRow++) {
			S.setElem(iRow, jCol, -S.getElem(iRow, jCol));
		}
	}

	Matrix ZKZt(K_);
	// make ZKZ'; it is a symmetric matrix
	ZKZt.premultZ(Z_);
	// will need ZK for calculating BLUPs; save it
	Matrix ZK(ZKZt);
	ZKZt.postmultZt(Z_);
	double offset = log(static_cast<double>(ZKZt.getNrows()));
	for (size_t ii = 0; ii < ZKZt.getNrows(); ii++) {
		double tmp = ZKZt.getElem(ii, ii) + offset;
		ZKZt.setElem(ii, ii, tmp);
	}

	S.symm('l', 'r', 1.0, ZKZt, 0.0, SHS);
	SHS.symm('l', 'r', 1.0, S, 0.0, u_); // u is now SHS

	vector<double> lam;
	S.resize(1, 1);
	u_.eigen('l', u_.getNcols() - XwIcpt.getNcols(), S, lam); // S now U(SHS)

	for (auto lamIt = lam.begin(); lamIt != lam.end(); ++lamIt) {
		if ((*lamIt) - offset <= 1e-10) {
			(*lamIt) = 0.0;
		} else {
			(*lamIt) -= offset;
		}

	}
	Matrix Eta;
	Y_.gemm(true, 1.0, S, false, 0.0, Eta); // U'Y

	for (size_t iRow = 0; iRow < Eta.getNrows(); iRow++) {
		for (size_t jCol = 0; jCol < Eta.getNcols(); jCol++) {
			Eta.setElem(iRow, jCol, pow2(Eta.getElem(iRow, jCol)));
		}
	}

	EmmREML reml(Eta, lam, 0);
	SHS = Y_; // replace SHS with Y. SHS will be modified in the loop
	for (size_t ii = 0; ii < ZKZt.getNrows(); ii++) {
		ZKZt.setElem(ii, ii, ZKZt.getElem(ii, ii) - offset);
	}
	vector<double> betaHat;
	Matrix xhx;  // X'H^{-1}X
	for (size_t jCol = 0; jCol < Y_.getNcols(); jCol++) {
		reml.setColID(jCol);
		double fMax = 0.0;
		maximizer(reml, 1e-6, delta_[jCol], fMax);

		u_ = ZKZt; // replacing Zu with ZKZ'
		for (size_t ii = 0; ii < u_.getNrows(); ii++) {
			u_.setElem(ii, ii, ZKZt.getElem(ii, ii) + delta_[jCol]);
		}
		u_.chol();
		u_.cholInv();          // now H^{-1}

		XwIcpt.symm('l', 'l', 1.0, u_, 0.0, S);        // S now H^{-1}X
		S.gemm(true, 1.0, XwIcpt, false, 0.0, xhx);    // xhx =  X'H^{-1}X
		xhx.chol();
		xhx.cholInv();                                 // xhx now (X'H^{-1}X)^{-1}
		S.gemc(true, 1.0, Y_, jCol, 0.0, betaHat);     // betaHat = X'H^{-1}Y[,jCol]
		beta_.setCol(jCol, betaHat);
		xhx.symc('l', 1.0, beta_, jCol, 0.0, betaHat);  // betaHat = (X'H^{-1}X)^{-1}X'H^{-1}Y[,jCol]
		beta_.setCol(jCol, betaHat);
		XwIcpt.gemc(false, 1.0, beta_, jCol, 0.0, betaHat); // betaHat now Xbeta; length nrow(X)

		for (size_t iRow = 0; iRow < Y_.getNrows(); iRow++) {
			betaHat[iRow] = Y_.getElem(iRow, jCol) - betaHat[iRow];  // betaHat = y - Xbeta
		}
		SHS.setCol(jCol, betaHat);
		betaHat.resize(0);                            // otherwise too long for gemc()
		u_.gemc(false, 1.0, SHS, jCol, 0.0, betaHat); // H^{-1}(y - Xbeta)
		SHS.setCol(jCol, betaHat);
		betaHat.resize(0);                            // essential for proper re-use at the top of the loop
	}
	SHS.gemm(true, 1.0, ZK, false, 0.0, u_); // replacing Zu with matrix of u

}

void MixedModel::hSq(vector<double> &out) const {
	out.resize(0);
	for (auto &dlt : delta_) {
		out.push_back(1.0/(1.0 + dlt));
	}
}

void MixedModel::gwa(){
	// calculate within-line means of Y
	Y_.premultZt(Z_);
	vector<double> nLn; // number of replicates per line
	Z_.colSums(nLn);
	Y_.colDivide(nLn);

	Matrix XwIcpt(1.0, Y_.getNrows(), 1);
	XwIcpt.appendCol(X_);
	Matrix uXb = u_;
	beta_.gemm(false, 1.0, XwIcpt, false, 1.0, uXb); // u = u + Xbeta
	Y_ = Y_ - uXb;

	// free memory from variables I no longer need
	nLn.resize(0);
	XwIcpt.resize(0, 0);
	uXb.resize(0,0);

	// determine the number of theards to use. There typically are two threads per core. It is worth doing no fewer than 200 SNPs per process.
	const size_t   Nsnp          = snps_->size()/K_.getNrows();
	const uint64_t nCores        = thread::hardware_concurrency()/2;
	const uint64_t minBlockSize  = 200;
	const uint64_t maxNumThreads = Nsnp/minBlockSize;
	uint64_t nSNPsDone           = 0;

	lPval_->resize(Y_.getNcols()*Nsnp);

	if ((maxNumThreads <= 1) || (nCores == 1)) {
		for (size_t iSNP = 0; iSNP < Nsnp; iSNP++) {
			oneSNP_(iSNP, Nsnp);
		}
	} else{
		uint64_t nThr      = (maxNumThreads <= nCores ? maxNumThreads : nCores);
		uint64_t blockSize = Nsnp/nThr;
		vector<thread> threads(nThr - 1);     // the main thread will do one block
		for (auto &thr : threads) {
			thr = thread{SNPblock(*this, nSNPsDone, blockSize)}; // has to be done on the fly because I deleted the copy constructor
		}
		if (nSNPsDone < Nsnp) { // if there are remaining SNPs to do, process the rest
			SNPblock block(*this, nSNPsDone, Nsnp-nSNPsDone);
			block();
		}

		for (auto &thr : threads) {
			if (thr.joinable()) {
				thr.join();
			}
		}
	}
}

void MixedModel::gwa(const uint32_t &nPer, vector<double> &fdr){
	this->gwa();
	const size_t   Nsnp          = snps_->size()/K_.getNrows();
	const uint64_t nCores        = thread::hardware_concurrency()/2;
	const uint64_t minBlockSize  = 200;
	const uint64_t maxNumThreads = Nsnp/minBlockSize;
	const uint64_t perOff        = Nsnp*nPer;

	vector<double>plPval(lPval_->size()*nPer);
	Matrix Ytmp(Y_);

	for (uint32_t i = 0; i < nPer; i++) {
		uint64_t nSNPsDone = 0;
		Y_ = Ytmp.rowShuffle();
		const uint64_t snpOff = Nsnp*i;
		if ((maxNumThreads <= 1) || (nCores == 1)) {
			for (size_t iSNP = 0; iSNP < Nsnp; iSNP++) {
				oneSNP_(iSNP, perOff, snpOff, plPval);
			}
		} else {
			uint64_t nThr      = (maxNumThreads <= nCores ? maxNumThreads : nCores);
			uint64_t blockSize = Nsnp/nThr;
			vector<thread> threads(nThr - 1);     // the main thread will do one block
			for (auto &thr : threads) {
				thr = thread{SNPblock(*this, nSNPsDone, blockSize, plPval, perOff, snpOff)}; // has to be done on the fly because I deleted the copy constructor
			}
			if (nSNPsDone < Nsnp) { // if there are remaining SNPs to do, process the rest
				SNPblock block(*this, nSNPsDone, Nsnp-nSNPsDone, plPval, perOff, snpOff);
				block();
			}

			for (auto &thr : threads) {
				if (thr.joinable()) {
					thr.join();
				}
			}
		}
	}
	fdr.resize(lPval_->size(), 1.0);
	vector<size_t> perOrder(plPval.size());
	vector<size_t> order(lPval_->size());

	// sort the real and permuted -lg(p) and calculate FDR for each trait
	for (size_t jTrait = 0; jTrait < Y_.getNcols(); jTrait++) {
		const size_t realStart = jTrait*Nsnp;
		const size_t realEnd   = realStart + Nsnp;
		const size_t perStart  = jTrait*perOff;
		const size_t perEnd    = perStart + perOff;
		quickSort(plPval,  perStart,  perEnd, perOrder);
		quickSort(*lPval_, realStart, realEnd, order);

		size_t curTopPer = perEnd - 1; // stores the index of the current top permuted -lg(p) that is smaller than the previous real data -lg(p); moving from the end of the perOrder vector
		double perCt     = 0.0;        // the number of permuted -lg(p)s larger than the current top real -lg(p)
		double realCt    = 1.0;
		// examine each real data -lg(p), starting with the largest (index is at the end of the order vector). Note that the trait offsets are already accounted for in order and perOrder vectors
		for (size_t iOrder = (realEnd - 1); iOrder >= realStart; iOrder--) {
			if (curTopPer == perStart) {    // reached the bottom of the permuted sample (in reverse order); any remaining FDR values are 1.0 (set at initialization).
				break;
			}
			double fVal  = 0.0;    // the F of Storey and Tibshirani (2003), equation [1] -- the number of false positives
			double sVal  = 0.0;    // the S of Storey and Tibshirani (2003), equation [1]

			const double curRealScore = (*lPval_)[ order[iOrder] ];

			// For the current real -lg(p), examine the sorted permuted -lg(p)s in order of decrease until we hit one that is smaller than the current real.
			// The number of permuted -lg(p)s larger than this one is the false positive count (all permuted assumed false)
			// Do not have to start from the beginning every time. All the permuted values tested with the previous real -lg(p) will be no smaller than the current one.
			double curPerScore = plPval[ perOrder[curTopPer] ]; // current largest permuted -lg(p)
			while (curRealScore <= curPerScore){
				curTopPer--;
				perCt += 1.0;
				curPerScore = plPval[ perOrder[curTopPer] ];
				if (curTopPer == perStart) {
					break;
				}
			}
			if (perCt > 0.0) {
				fVal = perCt/static_cast<double>(nPer); // scale by the number of permutations
				sVal = realCt;
				fdr[ order[iOrder] ] = (fVal > sVal ? 1.0 : fVal/sVal); // values > 1.0 can occur by chance but are not meaningful
			} else { // did not find any permuted -lg(p) larger than current real
				fdr[ order[iOrder] ] = 1.0/static_cast<double>(perOff);
			}
			realCt += 1.0;
			// This test is super importtant: since iOrder is unsigned, but I need the 0-th element, once I hit 0
			// the next decrement will wrap iOrder around, the test will be TRUE and we will continue with a crazy
			// out-of-bound value of iOrder, resulting in a segfault
			if (iOrder == 0) {
				break;
			}
		}
	}
}

void MixedModel::oneSNP_(const size_t &idx, const size_t &Nsnp){
	const size_t d     = Y_.getNcols();
	const size_t Nln   = Y_.getNrows();
	const size_t pad   = Nln*idx;      // where in the SNP vector the current SNP starts
	Matrix XtX(0.0, 2, 2);             // only lower triangle will be filled
	Matrix XtY(0.0, 2, d);             // for genotype == 0 || 1
	Matrix XtY2(0.0, 2, d);            // for genotype == 2
	// populate the X'X and X'Y matrices
	for (size_t iLn  = 0; iLn < Nln; iLn++) {
		const int32_t genotype = (*snps_)[iLn + pad];
		if (genotype == misTok_) {
			continue;
		}
		XtX.setElem(0, 0, XtX.getElem(0, 0) + 1.0);
		if (genotype == 0) {
			for (size_t jTrt = 0; jTrt < d; jTrt++) {
				XtY.setElem(0, jTrt, XtY.getElem(0, jTrt) + Y_.getElem(iLn, jTrt));
			}
		} else if (genotype == 1){
			XtX.setElem(1, 0, XtX.getElem(1, 0) + 1.0);
			XtX.setElem(1, 1, XtX.getElem(1, 1) + 1.0);
			for (size_t jTrt = 0; jTrt < d; jTrt++) {
				XtY.setElem(0, jTrt, XtY.getElem(0, jTrt) + Y_.getElem(iLn, jTrt));
				XtY.setElem(1, jTrt, XtY.getElem(1, jTrt) + Y_.getElem(iLn, jTrt));
			}
		} else if (genotype == 2){
			XtX.setElem(1, 0, XtX.getElem(1, 0) + 2.0);
			XtX.setElem(1, 1, XtX.getElem(1, 1) + 4.0);
			for (size_t jTrt = 0; jTrt < d; jTrt++) {
				XtY2.setElem(0, jTrt, XtY2.getElem(0, jTrt) + Y_.getElem(iLn, jTrt));
				XtY2.setElem(1, jTrt, XtY2.getElem(1, jTrt) + Y_.getElem(iLn, jTrt));
			}
		} else {
			throw("ERROR: unknown genotype in the SNP table in snpReg()");
		}
	}
	XtY2.rowMultiply(2.0, 1);
	XtY = XtY + XtY2;

	const size_t Npres = static_cast<size_t>(XtX.getElem(0, 0)); // this is the number of non-missing genotypes; save for future use
	const double df    = XtX.getElem(0, 0) - 2.0;                // degrees of freedom (n - q)
	XtX.chol();
	XtX.cholInv();
	XtY.symm('l', 'l', 1.0, XtX, 0.0, XtY2); // XtY2 is now the coefficient (beta) matrix
	XtY.resize(3, d); // there are three possible values of Xbeta, one for each genotype; store them in this matrix now
	for (size_t jTrt = 0; jTrt < d; jTrt++) {
		XtY.setElem(0, jTrt, XtY2.getElem(0, jTrt));
		XtY.setElem(1, jTrt, XtY2.getElem(0, jTrt) + XtY2.getElem(1, jTrt));
		XtY.setElem(2, jTrt, XtY2.getElem(0, jTrt) + 2.0*XtY2.getElem(1, jTrt));
	}
	Matrix rsdSq(Npres, d);
	size_t iPres = 0;
	for (size_t iLn = 0; iLn < Nln; iLn++) {
		const int genotype = (*snps_)[iLn + pad];
		if (genotype == misTok_) {
			continue;
		}
		for (size_t jTrt = 0; jTrt < d; jTrt++) {
			const double diff = Y_.getElem(iLn, jTrt) - XtY.getElem(genotype, jTrt);
			rsdSq.setElem(iPres, jTrt, pow2(diff));
		}
		iPres++;
	}
	vector<double> sSq;
	rsdSq.colSums(sSq);
	for (size_t jTrt = 0; jTrt < d; jTrt++) {
		const double sigSq = sSq[jTrt]*XtX.getElem(1, 1)/df;
		double fStat = pow2(XtY2.getElem(1, jTrt))/sigSq;
		fStat = df/(df + fStat);
		(*lPval_)[jTrt*Nsnp + idx] = -log10(betai(fStat, df/2.0, 0.5));
	}
}

void MixedModel::oneSNP_(const size_t &snpIdx, const size_t &perOff, const size_t &snpOff, vector<double> &lPval){
	const size_t d      = Y_.getNcols();
	const size_t Nln    = Y_.getNrows();
	const size_t pad    = Nln*snpIdx;        // position of the current SNP in the SNP vector (within a permutation)
	const size_t curSNP = snpOff + snpIdx;   // where in the current trait the SNP is, given the current permutation index
	Matrix XtX(0.0, 2, 2);                   // only lower triangle will be filled
	Matrix XtY(0.0, 2, d);                   // for genotype == 0 || 1
	Matrix XtY2(0.0, 2, d);                  // for genotype == 2
	// populate the X'X and X'Y matrices
	for (size_t iLn  = 0; iLn < Nln; iLn++) {
		const int32_t genotype = (*snps_)[iLn + pad];
		if (genotype == misTok_) {
			continue;
		}
		XtX.setElem(0, 0, XtX.getElem(0, 0) + 1.0);
		if (genotype == 0) {
			for (size_t jTrt = 0; jTrt < d; jTrt++) {
				XtY.setElem(0, jTrt, XtY.getElem(0, jTrt) + Y_.getElem(iLn, jTrt));
			}
		} else if (genotype == 1){
			XtX.setElem(1, 0, XtX.getElem(1, 0) + 1.0);
			XtX.setElem(1, 1, XtX.getElem(1, 1) + 1.0);
			for (size_t jTrt = 0; jTrt < d; jTrt++) {
				XtY.setElem(0, jTrt, XtY.getElem(0, jTrt) + Y_.getElem(iLn, jTrt));
				XtY.setElem(1, jTrt, XtY.getElem(1, jTrt) + Y_.getElem(iLn, jTrt));
			}
		} else if (genotype == 2){
			XtX.setElem(1, 0, XtX.getElem(1, 0) + 2.0);
			XtX.setElem(1, 1, XtX.getElem(1, 1) + 4.0);
			for (size_t jTrt = 0; jTrt < d; jTrt++) {
				XtY2.setElem(0, jTrt, XtY2.getElem(0, jTrt) + Y_.getElem(iLn, jTrt));
				XtY2.setElem(1, jTrt, XtY2.getElem(1, jTrt) + Y_.getElem(iLn, jTrt));
			}
		} else {
			throw("ERROR: unknown genotype in the SNP table in snpReg()");
		}
	}
	XtY2.rowMultiply(2.0, 1);
	XtY = XtY + XtY2;

	const size_t Npres = static_cast<size_t>(XtX.getElem(0, 0)); // this is the number of non-missing genotypes; save for future use
	const double df    = XtX.getElem(0, 0) - 2.0;                // degrees of freedom (n - q)
	XtX.chol();
	XtX.cholInv();
	XtY.symm('l', 'l', 1.0, XtX, 0.0, XtY2); // XtY2 is now the coefficient (beta) matrix
	XtY.resize(3, d); // there are three possible values of Xbeta, one for each genotype; store them in this matrix now
	for (size_t jTrt = 0; jTrt < d; jTrt++) {
		XtY.setElem(0, jTrt, XtY2.getElem(0, jTrt));
		XtY.setElem(1, jTrt, XtY2.getElem(0, jTrt) + XtY2.getElem(1, jTrt));
		XtY.setElem(2, jTrt, XtY2.getElem(0, jTrt) + 2.0*XtY2.getElem(1, jTrt));
	}
	Matrix rsdSq(Npres, d);
	size_t iPres = 0;
	for (size_t iLn = 0; iLn < Nln; iLn++) {
		const int genotype = (*snps_)[iLn + pad];
		if (genotype == misTok_) {
			continue;
		}
		for (size_t jTrt = 0; jTrt < d; jTrt++) {
			const double diff = Y_.getElem(iLn, jTrt) - XtY.getElem(genotype, jTrt);
			rsdSq.setElem(iPres, jTrt, pow2(diff));
		}
		iPres++;
	}
	vector<double> sSq;
	rsdSq.colSums(sSq);
	for (size_t jTrt = 0; jTrt < d; jTrt++) {
		const double sigSq = sSq[jTrt]*XtX.getElem(1, 1)/df;
		double fStat = pow2(XtY2.getElem(1, jTrt))/sigSq;
		fStat = df/(df + fStat);
		lPval[perOff*jTrt + curSNP] = -log10(betai(fStat, df/2.0, 0.5));
	}
}

// SNPblock methods
void SNPblock::operator()(){
	if (plPval_ == nullptr) {
		const size_t Nsnp = mmObj_->snps_->size()/(mmObj_->Y_.getNrows());
		for (size_t kSNP = blockStart_; kSNP < blockStart_ + blockSize_; kSNP++) {
			mmObj_->oneSNP_(kSNP, Nsnp);
		}
	} else {
		for (size_t kSNP = blockStart_; kSNP < blockStart_ + blockSize_; kSNP++) {
			mmObj_->oneSNP_(kSNP, perOff_, snpOff_, *plPval_);
		}
	}
}


