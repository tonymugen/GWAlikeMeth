//
//  likeMeth.cpp
//  QuaGen
//
//  Created by Tony Mugen on 3/29/17.
//

/// Likelihood methods for quantitative genetics
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2017 Anthony J. Greenberg
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
	for (size_t s = 0; s < _etaSq->getNrows(); s++) {
		ld    = (*_lambda)[s] + delta;
		frac += _etaSq->getElem(s, _jCol)/ld;
		llam += log(ld);
	}
	// final likelihood
	double llik = -static_cast<double>(_lambda->size())*log(frac) - llam;
	return llik;

}

void SNPblock::operator()(){
	const size_t Nrow = _rsp->getNrows();
	const size_t Ncol = _rsp->getNcols();


	for (size_t kSNP = _blockStart; kSNP < _blockStart + _blockSize; kSNP++) {
		snpReg(*_rsp, _snp + kSNP*Nrow, -9, _lPval + kSNP*Ncol);
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


