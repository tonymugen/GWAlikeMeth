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

#include "likeMeth.hpp"
#include "utilities.hpp"

using std::cerr;
using std::endl;
using std::vector;
using locMatrix::Matrix;

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

void SNPblock::operator()(){
	const size_t Nrow = rsp_->getNrows();
	const size_t Ncol = rsp_->getNcols();


	for (size_t kSNP = blockStart_; kSNP < blockStart_ + blockSize_; kSNP++) {
		snpReg(*rsp_, snp_ + kSNP*Nrow, -9, lPval_ + kSNP*Ncol);
	}
}

void solveMM(const Matrix &Y, const Matrix &K, const Matrix &Z, Matrix &u, vector<double> &mu){
	Matrix ZKZt(K);
	if (mu.size()) {
		mu.resize(0);
	}
	double invN = -1.0/static_cast<double>(Z.getNrows());
	Matrix S(invN, Z.getNrows(), Z.getNrows());
	for (size_t i = 0; i < Z.getNrows(); i++) {
		S.setElem(i, i, 1.0 + invN);
	}

	// make ZKZ'; it is a symmetric matrix
	ZKZt.premultZ(Z);
	// will need ZK for calculating BLUPs; save it
	Matrix ZK(ZKZt);
	ZKZt.postmultZt(Z);
	double offset = log(static_cast<double>(ZKZt.getNrows()));
	for (size_t ii = 0; ii < ZKZt.getNrows(); ii++) {
		double tmp = ZKZt.getElem(ii, ii) + offset;
		ZKZt.setElem(ii, ii, tmp);
	}

	Matrix SHS; // SHS matrix (equation 5 of Kang et al.)
	S.symm('l', 'r', 1.0, ZKZt, 0.0, SHS);
	SHS.symm('l', 'r', 1.0, S, 0.0, u); // u is now SHS

	vector<double> lam;
	S.resize(1, 1);
	u.eigen('l', u.getNcols() - 1, S, lam); // S now U(SHS)

	for (auto lamIt = lam.begin(); lamIt != lam.end(); ++lamIt) {
		if ((*lamIt) - offset <= 1e-10) {
			(*lamIt) = 0.0;
		} else {
			(*lamIt) -= offset;
		}

	}

	Matrix Eta;
	Y.gemm(true, 1.0, S, false, 0.0, Eta); // U'Y

	for (size_t iRow = 0; iRow < Eta.getNrows(); iRow++) {
		for (size_t jCol = 0; jCol < Eta.getNcols(); jCol++) {
			Eta.setElem(iRow, jCol, pow2(Eta.getElem(iRow, jCol)));
		}
	}

	EmmREML reml(Eta, lam, 0);
	SHS = Y; // replace SHS with Y. SHS will be modified in the loop
	vector<double> XtH(ZKZt.getNrows());
	for (size_t ii = 0; ii < ZKZt.getNrows(); ii++) {
		ZKZt.setElem(ii, ii, ZKZt.getElem(ii, ii) - offset);
	}

	for (size_t jCol = 0; jCol < Y.getNcols(); jCol++) {
		reml.setColID(jCol);
		double deltaOpt = 0.0;
		double fMax     = 0.0;
		maximizer(reml, 1e-6, deltaOpt, fMax);

		u = ZKZt; // replacing Zu with ZKZ'
		for (size_t ii = 0; ii < u.getNrows(); ii++) {
			u.setElem(ii, ii, ZKZt.getElem(ii, ii) + deltaOpt);
		}
		u.chol();
		u.cholInv();          // now H^{-1}
		u.colSums(XtH);       // X'H^{-1}, since X is a column of 1s (intercept)
		double XtHX    = 0.0; // X'H^{-1}X
		double betaLoc = 0.0;
		for (size_t iRow = 0; iRow < Y.getNrows(); iRow++) {
			XtHX    += XtH[iRow];
			betaLoc += Y.getElem(iRow, jCol) * XtH[iRow];
		}

		betaLoc /= XtHX;
		mu.push_back(betaLoc);

		SHS.colSub(betaLoc, jCol);
		u.gemc(false, 1.0, SHS, jCol, 0.0, XtH);
		SHS.setCol(jCol, XtH);

	}
	SHS.gemm(true, 1.0, ZK, false, 0.0, u); // replacing Zu with matrix of u

}
void solveMM(const Matrix &Y, const Matrix &K, const Matrix &Z, const Matrix &X, Matrix &u, Matrix &beta){
	if (X.getNcols() >= X.getNrows()) {
		cerr << "ERROR: covariate matrix X (" << X.getNrows() << " x " << X.getNcols() << ") is not full rank" << endl;
		return;
	}
	if ( (beta.getNcols() != Y.getNcols()) || (beta.getNrows() != (X.getNcols() + 1)) ) {
		beta.resize(X.getNcols() + 1, Y.getNcols());
	}
	// add the intercept column to X
	Matrix XwIcpt(1.0, X.getNrows(), 1);
	XwIcpt.appendCol(X);
	if (XwIcpt.getNrows() == K.getNrows()) {
		XwIcpt.premultZ(Z);
	} else if (XwIcpt.getNrows() != Y.getNrows()) {
		cerr << "ERROR: number of rows in X (" << X.getNrows() << ") not equal to number of lines (" << K.getNrows() << ") or number of data points (" << Y.getNrows() << ")" << endl;
		return;
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

	Matrix ZKZt(K);
	// make ZKZ'; it is a symmetric matrix
	ZKZt.premultZ(Z);
	// will need ZK for calculating BLUPs; save it
	Matrix ZK(ZKZt);
	ZKZt.postmultZt(Z);
	double offset = log(static_cast<double>(ZKZt.getNrows()));
	for (size_t ii = 0; ii < ZKZt.getNrows(); ii++) {
		double tmp = ZKZt.getElem(ii, ii) + offset;
		ZKZt.setElem(ii, ii, tmp);
	}

	S.symm('l', 'r', 1.0, ZKZt, 0.0, SHS);
	SHS.symm('l', 'r', 1.0, S, 0.0, u); // u is now SHS

	vector<double> lam;
	S.resize(1, 1);
	u.eigen('l', u.getNcols() - XwIcpt.getNcols(), S, lam); // S now U(SHS)

	for (auto lamIt = lam.begin(); lamIt != lam.end(); ++lamIt) {
		if ((*lamIt) - offset <= 1e-10) {
			(*lamIt) = 0.0;
		} else {
			(*lamIt) -= offset;
		}

	}
	Matrix Eta;
	Y.gemm(true, 1.0, S, false, 0.0, Eta); // U'Y

	for (size_t iRow = 0; iRow < Eta.getNrows(); iRow++) {
		for (size_t jCol = 0; jCol < Eta.getNcols(); jCol++) {
			Eta.setElem(iRow, jCol, pow2(Eta.getElem(iRow, jCol)));
		}
	}

	EmmREML reml(Eta, lam, 0);
	SHS = Y; // replace SHS with Y. SHS will be modified in the loop
	for (size_t ii = 0; ii < ZKZt.getNrows(); ii++) {
		ZKZt.setElem(ii, ii, ZKZt.getElem(ii, ii) - offset);
	}
	vector<double> betaHat;
	Matrix xhx;  // X'H^{-1}X
	for (size_t jCol = 0; jCol < Y.getNcols(); jCol++) {
		reml.setColID(jCol);
		double deltaOpt = 0.0;
		double fMax     = 0.0;
		maximizer(reml, 1e-6, deltaOpt, fMax);

		u = ZKZt; // replacing Zu with ZKZ'
		for (size_t ii = 0; ii < u.getNrows(); ii++) {
			u.setElem(ii, ii, ZKZt.getElem(ii, ii) + deltaOpt);
		}
		u.chol();
		u.cholInv();          // now H^{-1}

		XwIcpt.symm('l', 'l', 1.0, u, 0.0, S);         // S now H^{-1}X
		S.gemm(true, 1.0, XwIcpt, false, 0.0, xhx);    // xhx =  X'H^{-1}X
		xhx.chol();
		xhx.cholInv();                                 // xhx now (X'H^{-1}X)^{-1}
		S.gemc(true, 1.0, Y, jCol, 0.0, betaHat);      // betaHat = X'H^{-1}Y[,jCol]
		beta.setCol(jCol, betaHat);
		xhx.symc('l', 1.0, beta, jCol, 0.0, betaHat);  // betaHat = (X'H^{-1}X)^{-1}X'H^{-1}Y[,jCol]
		beta.setCol(jCol, betaHat);
		XwIcpt.gemc(false, 1.0, beta, jCol, 0.0, betaHat); // betaHat now Xbeta; length nrow(X)

		for (size_t iRow = 0; iRow < Y.getNrows(); iRow++) {
			betaHat[iRow] = Y.getElem(iRow, jCol) - betaHat[iRow];  // betaHat = y - Xbeta
		}
		SHS.setCol(jCol, betaHat);
		betaHat.resize(0);                           // otherwise too long for gemc()
		u.gemc(false, 1.0, SHS, jCol, 0.0, betaHat); // H^{-1}(y - Xbeta)
		SHS.setCol(jCol, betaHat);
		betaHat.resize(0);                           // essential for proper re-use at the top of the loop
	}
	SHS.gemm(true, 1.0, ZK, false, 0.0, u); // replacing Zu with matrix of u

}

void snpReg(const Matrix &Y, const int *snp, const int &misTok, double *lPval){
	const size_t d   = Y.getNcols();
	const size_t Nln = Y.getNrows();
	Matrix XtX(0.0, 2, 2);  // only lower triangle will be filled
	Matrix XtY(0.0, 2, d);  // for genotype == 0 || 1
	Matrix XtY2(0.0, 2, d); // for genotype == 2
	// populate the X'X and X'Y matrices
	for (size_t iLn  = 0; iLn < Nln; iLn++) {
		const int genotype = snp[iLn];
		if (genotype == misTok) {
			continue;
		}
		XtX.setElem(0, 0, XtX.getElem(0, 0) + 1.0);
		if (genotype == 0) {
			for (size_t jTrt = 0; jTrt < d; jTrt++) {
				XtY.setElem(0, jTrt, XtY.getElem(0, jTrt) + Y.getElem(iLn, jTrt));
			}
		} else if (genotype == 1){
			XtX.setElem(1, 0, XtX.getElem(1, 0) + 1.0);
			XtX.setElem(1, 1, XtX.getElem(1, 1) + 1.0);
			for (size_t jTrt = 0; jTrt < d; jTrt++) {
				XtY.setElem(0, jTrt, XtY.getElem(0, jTrt) + Y.getElem(iLn, jTrt));
				XtY.setElem(1, jTrt, XtY.getElem(1, jTrt) + Y.getElem(iLn, jTrt));
			}
		} else if (genotype == 2){
			XtX.setElem(1, 0, XtX.getElem(1, 0) + 2.0);
			XtX.setElem(1, 1, XtX.getElem(1, 1) + 4.0);
			for (size_t jTrt = 0; jTrt < d; jTrt++) {
				XtY2.setElem(0, jTrt, XtY2.getElem(0, jTrt) + Y.getElem(iLn, jTrt));
				XtY2.setElem(1, jTrt, XtY2.getElem(1, jTrt) + Y.getElem(iLn, jTrt));
			}
		} else {
			cerr << "ERROR: unknown genotype " << genotype << " in the SNP table in snpReg()" << endl;
			exit(100);
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
		const int genotype = snp[iLn];
		if (genotype == misTok) {
			continue;
		}
		for (size_t jTrt = 0; jTrt < d; jTrt++) {
			const double diff = Y.getElem(iLn, jTrt) - XtY.getElem(genotype, jTrt);
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
		lPval[jTrt] = -log10(betai(fStat, df/2.0, 0.5));
	}

}

// MixedModel methods

MixedModel::MixedModel(const vector<double> &yvec, const vector<double> &kvec, const vector<double> &repFac, const size_t &d, const size_t &Ngen, const size_t &N) : Y_{Matrix(yvec, N, d)}, K_{Matrix(kvec, Ngen, Ngen)}, Z_{Matrix(repFac, N, Ngen)} {
	Matrix ZKZt(K_);
	vector<double> mu;
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

	for (size_t jCol = 0; jCol < Y_.getNcols(); jCol++) {
		reml.setColID(jCol);
		double deltaOpt = 0.0;
		double fMax     = 0.0;
		maximizer(reml, 1e-6, deltaOpt, fMax);

		u_ = ZKZt; // replacing Zu with ZKZ'
		for (size_t ii = 0; ii < u_.getNrows(); ii++) {
			u_.setElem(ii, ii, ZKZt.getElem(ii, ii) + deltaOpt);
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
		mu.push_back(betaLoc);

		SHS.colSub(betaLoc, jCol);
		u_.gemc(false, 1.0, SHS, jCol, 0.0, XtH);
		SHS.setCol(jCol, XtH);

	}
	SHS.gemm(true, 1.0, ZK, false, 0.0, u_); // replacing Zu with matrix of u


}

MixedModel::MixedModel(const vector<double> &yvec, const vector<double> &kvec, const vector<double> &repFac, const vector<double> &xvec, const size_t &d, const size_t &Ngen, const size_t &N) : Y_{Matrix(yvec, N, d)}, K_{Matrix(kvec, Ngen, Ngen)}, Z_{Matrix(repFac, N, Ngen)}, X_{Matrix(xvec, N, xvec.size()/N)} {

	if (X_.getNcols() >= X_.getNrows()) {
		throw string("ERROR: covariate matrix X_ in MixedModel constructor is not full rank");
	}
	if ( (beta_.getNcols() != Y_.getNcols()) || (beta_.getNrows() != (X_.getNcols() + 1)) ) {
		beta_.resize(X_.getNcols() + 1, Y_.getNcols());
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
		double deltaOpt = 0.0;
		double fMax     = 0.0;
		maximizer(reml, 1e-6, deltaOpt, fMax);

		u_ = ZKZt; // replacing Zu with ZKZ'
		for (size_t ii = 0; ii < u_.getNrows(); ii++) {
			u_.setElem(ii, ii, ZKZt.getElem(ii, ii) + deltaOpt);
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


