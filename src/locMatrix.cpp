/*
 * Copyright (c) 2016 Anthony J. Greenberg
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


/// C++ matrix class for development.
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2016 Anthony J. Greenberg
 * \version 0.1
 *
 * This is the class implementation file for the experimental Matrix class. This version is for including in R packages, so it uses the R BLAS and LAPACK interfaces.
 *
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <climits>
#include <utility>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

#include "locMatrix.hpp"
#include "random.hpp"

using std::endl;
using std::flush;
using std::ofstream;
using std::ifstream;
using std::stringstream;
using std::ios;
using std::vector;
using std::string;
using std::stod;
using std::fill;
using std::memcpy;
using std::nan;
using std::numeric_limits;

using namespace BayesicSpace;


// Friend functions
Matrix BayesicSpace::operator*(const double &scal, const Matrix &m){
	Matrix res(m);
	for (size_t iElm = 0; iElm < m.Ncol_*m.Nrow_; iElm++) {
		res.data_[iElm] *= scal;
	}

	return res;
}

Matrix BayesicSpace::operator+(const double &scal, const Matrix &m){
	Matrix res(m);
	for (size_t iElm = 0; iElm < m.Ncol_*m.Nrow_; iElm++) {
		res.data_[iElm] += scal;
	}

	return res;
}

// Matrix methods

Matrix::Matrix(const size_t &nrow, const size_t &ncol) : Nrow_(nrow), Ncol_(ncol) { // not using uniform initialization list because Xcode code completion goes nuts
	if ( (Nrow_ > INT_MAX) || (Ncol_ > INT_MAX) ) {
		throw string("ERROR: Matrix dimensions exceed INT_MAX in Matrix initializing constructor");
	}
	size_t Ntot = Nrow_ * Ncol_;
	if (!Ntot) { // one of the dimensions is zero
		throw string("ERROR: one of the dimensions is 0 in the Matrix initializing constructor");
	}
	if ((Ntot < Nrow_) || (Ntot < Ncol_)) { // this happens only if there is wrap-around
		throw string("ERROR: dimensions are too large resulting in wrap-around in the Matrix initializing constructor");
	}
	data_ = new double[Ntot];

}

// Using the C++-11 delegating constructor facility
Matrix::Matrix(const double &val, const size_t &nrow, const size_t &ncol) : Matrix(nrow, ncol) {
	size_t Ntot = Nrow_ * Ncol_;
	fill(data_, data_ + Ntot, val);

}

Matrix::Matrix(const double inArr[], const size_t &nrow, const size_t &ncol) : Matrix(nrow, ncol) {
	size_t Ntot = Nrow_ * Ncol_;
	memcpy(data_, inArr, Ntot*sizeof(*inArr));
}

Matrix::Matrix(const vector<double> &inVec, const size_t &nrow, const size_t &ncol) : Matrix(nrow, ncol) {
	size_t Ntot = Nrow_ * Ncol_;

	if (inVec.size() < Ntot){
		throw string("ERROR: input vector is not long enough in the Matrix vector-based constructor");
	}

	memcpy(data_, inVec.data(), Ntot*sizeof(double));
}

Matrix::Matrix(const string &fileName, const char &delim): Nrow_(0), Ncol_(0) {
	ifstream matIn(fileName.c_str());

	if (!matIn) {
		string errMsg = "ERROR: cannot open file " + fileName + " in Matrix constructor from file";
		throw errMsg;
	}

	// Copy all the lines into a vector of strings.
	// Cannot assign the values to the data_ array on the fly because I need to know the dimensions to allocate the data_ array.
	vector<string> lines;
	string currentElem;

	while (getline(matIn, currentElem)){
		lines.push_back(currentElem);
	}
	matIn.close();

	Nrow_ = lines.size();
	if (Nrow_ == 0) {
		string errMsg = "ERROR: no rows in the file " + fileName + "ERROR: no rows in the file";
		throw errMsg;
	}

	// going through the first line is enough to establish
	vector<double> firstLine;
	stringstream lineStream;
	lineStream.str(lines[0]);
	while (getline(lineStream, currentElem, delim)) {
		firstLine.push_back(stod(currentElem));
	}
	lineStream.str("");

	Ncol_ = firstLine.size();
	if (Ncol_ == 0) {
		string errMsg = "ERROR: no rows in the file " + fileName + " in Matrix constructor from file";
		throw errMsg;
	}

	size_t Ntot = Ncol_ * Nrow_;
	if ((Ntot < Nrow_) || (Ntot < Ncol_)) { // this happens only if there is wrap-around
		throw string("ERROR: dimensions are too large resulting in wrap-around in the Matrix constructor from file");
	}
	data_ = new double[Ntot];

	size_t jCol = 0;
	for (auto flIt = firstLine.begin(); flIt != firstLine.end(); ++flIt) {
		data_[Nrow_*jCol] = *flIt; // just the first row
		jCol++;
	}
	firstLine.resize(0);

	// now go through the rest of the lines
	size_t iRow = 1;
	for (auto lnIt = lines.begin() + 1; lnIt != lines.end(); ++lnIt) {
		stringstream elemStrm;
		elemStrm.str(*lnIt);
		jCol = 0;
		while (getline(elemStrm, currentElem, delim)) {
			data_[Nrow_*jCol + iRow] = stod(currentElem);
			jCol++;
		}
		iRow++;
		if (jCol != Ncol_) {
			string errMsg = "ERROR: in file " + fileName + " a not all lines have the same number of columns (Matrix constructor from file)";
			throw errMsg;
		}
	}
}

Matrix::~Matrix(){
	delete [] data_;
	data_ = nullptr;
	Ncol_ = 0;
	Nrow_ = 0;

}

Matrix::Matrix(const Matrix &inMat){
	Ncol_ = inMat.Ncol_;
	Nrow_ = inMat.Nrow_;
	if (Ncol_ && Nrow_) {
		data_ = new double[Nrow_ * Ncol_];
		memcpy(data_, inMat.data_, (Ncol_ * Nrow_)*sizeof(double));
	} else {
		data_ = nullptr;
	}
}

Matrix& Matrix::operator=(const Matrix &inMat){
	if (this != &inMat) {
		delete [] data_;
		Ncol_ = inMat.Ncol_;
		Nrow_ = inMat.Nrow_;
		if (Ncol_ && Nrow_) {
			data_ = new double[Nrow_ * Ncol_];
			memcpy(data_, inMat.data_, (Ncol_ * Nrow_)*sizeof(double));
		} else {
			data_ = nullptr;
		}
	}

	return *this;
}
Matrix::Matrix(Matrix &&inMat) {
	data_ = inMat.data_;
	Ncol_ = inMat.Ncol_;
	Nrow_ = inMat.Nrow_;

	inMat.data_ = nullptr;
	inMat.Ncol_ = 0;
	inMat.Nrow_ = 0;
}
Matrix& Matrix::operator=(Matrix &&inMat){
	if (this != &inMat) {
		delete [] data_;
		data_ = inMat.data_;
		Ncol_ = inMat.Ncol_;
		Nrow_ = inMat.Nrow_;

		inMat.data_ = nullptr;
		inMat.Ncol_ = 0;
		inMat.Nrow_ = 0;
	}

	return *this;
}

double Matrix::getElem(const size_t& iRow, const size_t &jCol) const{
#ifndef LMRG_CHECK_OFF
	if ((iRow >= Nrow_) || (jCol >= Ncol_)) {
		throw("ERROR: element out of range in getElem()");
	}
#endif

	return data_[Nrow_*jCol + iRow];
}

void Matrix::setElem(const size_t& iRow, const size_t &jCol, const double &input){
#ifndef LMRG_CHECK_OFF
	if ((iRow >= Nrow_) || (jCol >= Ncol_)) {
		throw string("ERROR: element out of range in setElem()");
	}
#endif

	data_[Nrow_*jCol + iRow] = input;

}

void Matrix::setCol(const size_t jCol, const vector<double> data){
#ifndef LMRG_CHECK_OFF
	if (jCol >= Ncol_) {
		throw string("ERROR: column index out of range in setCol()");
	}
	if (data.size() < Nrow_) {
		throw string("ERROR: vector length smaller than the number of rows in setCol()");
	}
#endif
	double *colBeg = data_ + jCol*Nrow_;
	memcpy(colBeg, data.data(), Nrow_*sizeof(double));

}

void Matrix::resize(const size_t &nrow, const size_t &ncol){
	delete [] data_;
	Nrow_ = nrow;
	Ncol_ = ncol;
	if (Ncol_ && Nrow_) {
		data_ = new double[Nrow_ * Ncol_]();
	} else {
		data_ = nullptr;
	}

}

void Matrix::save(const string &outFileName) const {
	remove(outFileName.c_str());
	ofstream outFl(outFileName.c_str(), ios::app);

	if (!outFl) {
		string errMsg = "ERROR: cannot open file " + outFileName + " for writing to save Matrix";
		throw errMsg;
	}
	if ((Nrow_ == 0) || (Ncol_ == 0)) {
		outFl.close();
		return;
	}
	for (size_t iRow = 0; iRow < Nrow_; iRow++) {
		// first element
		outFl << data_[iRow] << flush;
		for (size_t jCol = 1; jCol < Ncol_; jCol++) {
			outFl << "\t" << data_[Nrow_*jCol + iRow] << flush;
		}
		outFl << endl;
	}

	outFl.close();
}

void Matrix::vectorize(vector<double> &out) const {
	out.resize(Nrow_*Ncol_);

	memcpy(out.data(), data_, out.size()*sizeof(double));
}

void Matrix::chol(){
#ifndef LMRG_CHECK_OFF
	if (Nrow_ != Ncol_) {
		throw string("ERROR: matrix has to be symmetric for Cholesky decomposition");
	}
	if (Nrow_ > INT_MAX) {
		throw string("ERROR: matrix dimension too big to safely convert to int in in-place Cholesky decomposition");
	}
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw string("ERROR: one of the dimensions is zero");
	}
#endif
	int info = 0;
	char tri = 'L';

	int N = static_cast<int>(Nrow_); // conversion should be OK: magnitude of Nrow_ checked in the constructor
	dpotrf_(&tri, &N, data_, &N, &info);
	if (info < 0) {
		throw string("ERROR: illegal element in in-place Cholesky decomposition");
	} else if (info > 0) {
		throw string("ERROR: matrix is not positive definite in in-place Cholesky decomposition");
	}

}

void Matrix::chol(Matrix &out) const {
#ifndef LMRG_CHECK_OFF
	if (Nrow_ != Ncol_) {
		throw string("ERROR: matrix has to be symmetric for Cholesky decomposition");
	}
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw string("ERROR: one of the dimensions is zero");
	}
#endif

	if ((Nrow_ != out.Nrow_) || (Ncol_ != out.Ncol_)) {
		out.resize(Nrow_, Ncol_);
	}
	memcpy(out.data_, data_, (Nrow_ * Ncol_)*sizeof(double));

	int info = 0;
	char tri = 'L';

	int N = static_cast<int>(Nrow_); // conversion should be safe: Nrow_ magnitude checked during construction
	dpotrf_(&tri, &N, out.data_, &N, &info);
	if (info < 0) {
		throw string("ERROR: illegal matrix element in copy Cholesky decomposition");
	} else if (info > 0) {
		throw string("ERROR: matrix is not positive definite in copy Cholesky decomposition");
	}

}

void Matrix::cholInv(){
#ifndef LMRG_CHECK_OFF
	if (Nrow_ != Ncol_) {
		throw string("ERROR: matrix has to be symmetric for Cholesky inversion");
	}
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw string("ERROR: one of the dimensions is zero");
	}
#endif
	int info = 0;
	char tri = 'L';

	int N = static_cast<int>(Nrow_); // conversion should be safe: Nrow_ magnitude checked during construction
	dpotri_(&tri, &N, data_, &N, &info);
	if (info < 0) {
		throw string("ERROR: illegal matrix element in-place Cholesky inversion");
	} else if (info > 0) {
		throw string("ERROR: a diagonal element of the matrix is zero. Cannot complete in-place Cholesky inversion");
	}
	// copying the lower triangle to the upper
	for (size_t iRow = 0; iRow < Nrow_; iRow++) {
		for (size_t jCol = 0; jCol < iRow; jCol++) {
			data_[Nrow_*iRow + jCol] = data_[Nrow_*jCol + iRow];
		}
	}

}

void Matrix::cholInv(Matrix &out) const {
#ifndef LMRG_CHECK_OFF
	if (Nrow_ != Ncol_) {
		throw string("ERROR: matrix has to be square for Cholesky inversion");
	}
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw string("ERROR: one of the dimensions is zero");
	}
#endif

	if ((Nrow_ != out.Nrow_) || (Ncol_ != out.Ncol_)) {
		out.resize(Nrow_, Ncol_);
	}
	memcpy(out.data_, data_, (Nrow_ * Ncol_)*sizeof(double));

	int info = 0;
	char tri = 'L';

	int N = static_cast<int>(Nrow_); // safe to convert: Nrow_ checked at construction
	dpotri_(&tri, &N, out.data_, &N, &info);
	if (info < 0) {
		throw string("ERROR: illegal matrix element in copy Cholesky inversion");
	} else if (info > 0) {
		throw string("ERROR: a diagonal element of the matrix is zero. Cannot complete copy Cholesky inversion");
	}
	for (size_t iRow = 0; iRow < Nrow_; iRow++) {
		for (size_t jCol = 0; jCol < iRow; jCol++) {
			out.data_[Nrow_*iRow + jCol] = out.data_[Nrow_*jCol + iRow];
		}
	}
}

void Matrix::svd(Matrix &U, vector<double> &s){
#ifndef LMRG_CHECK_OFF
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw string("ERROR: one of the dimensions is zero");
	}
#endif

	if ((Nrow_ != U.Nrow_) || (U.Nrow_ != U.Ncol_)) {
		U.resize(Nrow_, Nrow_);
	}
	if (s.size() < Ncol_) {
		s.resize(Ncol_, 0.0);
	}
	int Nvt = 1;
	vector<double>vt(1, 0.0);
	int resSVD = 0;
	int Nw = -1;    // set this to pre-run dgesvd_ for calculation of workspace size
	vector<double>workArr(1, 0.0);
	char jobu  = 'A';
	char jobvt = 'N';
	// the following casts are safe because dimensions are checked at construction
	int Nr = static_cast<int>(Nrow_);
	int Nc = static_cast<int>(Ncol_);

	// first calculate working space
	dgesvd_(&jobu, &jobvt, &Nr, &Nc, data_, &Nr, s.data(), U.data_, &Nr, vt.data(), &Nvt, workArr.data(), &Nw, &resSVD);
	Nw = workArr[0];
	workArr.resize(Nw, 0.0);
	dgesvd_(&jobu, &jobvt, &Nr, &Nc, data_, &Nr, s.data(), U.data_, &Nr, vt.data(), &Nvt, workArr.data(), &Nw, &resSVD);
	workArr.resize(0);
	if (resSVD < 0) {
		throw string("ERROR: illegal matrix element in SVD");
	} else if (resSVD > 0){
		throw string("ERROR: DBDSQR did not converge in SVD");
	}

}

void Matrix::svdSafe(Matrix &U, vector<double> &s) const {
#ifndef LMRG_CHECK_OFF
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw string("ERROR: one of the dimensions is zero");
	}
#endif

	double *dataCopy = new double[Nrow_ * Ncol_];
	memcpy(dataCopy, data_, (Nrow_ * Ncol_)*sizeof(double));

	if ((Nrow_ != U.Nrow_) || (U.Nrow_ != U.Ncol_)) {
		U.resize(Nrow_, Nrow_);
	}
	if (s.size() < Ncol_) {
		s.resize(Ncol_, 0.0);
	}
	int Nvt = 1;
	vector<double>vt(1, 0.0);
	int resSVD = 0;
	int Nw = -1;    // set this to pre-run dgesvd_ for calculation of workspace size
	vector<double>workArr(1, 0.0);
	char jobu  = 'A';
	char jobvt = 'N';
	// the folloeing casts are safe because the dimensions are checked at construction
	int Nr = static_cast<int>(Nrow_);
	int Nc = static_cast<int>(Ncol_);

	// first calculate working space
	dgesvd_(&jobu, &jobvt, &Nr, &Nc, dataCopy, &Nr, s.data(), U.data_, &Nr, vt.data(), &Nvt, workArr.data(), &Nw, &resSVD);
	Nw = workArr[0];
	workArr.resize(Nw, 0.0);
	dgesvd_(&jobu, &jobvt, &Nr, &Nc, dataCopy, &Nr, s.data(), U.data_, &Nr, vt.data(), &Nvt, workArr.data(), &Nw, &resSVD);
	workArr.resize(0);
	if (resSVD < 0) {
		throw string("ERROR: illegal matrix element in safe SVD");
		exit(14);
	} else if (resSVD > 0){
		throw string("ERROR: DBDSQR did not converge in safe SVD");
	}
	delete [] dataCopy;
}

void Matrix::eigen(const char &tri, Matrix &U, vector<double> &lam){
#ifndef LMRG_CHECK_OFF
	if (Nrow_ != Ncol_) {
		throw string("ERROR: matrix has to be at least square in eigen()");
	}
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw string("ERROR: one of the dimensions is zero");
	}
#endif

	// test the output size and adjust if necessary
	if ((Ncol_ > U.Nrow_) || (Ncol_ > U.Ncol_)) {
		U.resize(Ncol_, Ncol_);
	}
	if (Ncol_ > lam.size()) {
		lam.resize(Ncol_, 0.0);
	}

	char jobz  = 'V'; // computing eigenvectors
	char range = 'A'; // doing all of them
	char uplo;
	if (tri == 'u') {
		uplo = 'U';
	} else if (tri == 'l'){
		uplo = 'L';
	} else {
		throw string("ERROR: unknown triangle indicator in eigen()");
	}
	// the following casts are safe because Nrow_ magnitude is checked at construction
	int N   = static_cast<int>(Nrow_);
	int lda = static_cast<int>(Nrow_);
	// placeholder variables. Not referenced since we are computing all eigenvectors
	double vl = 0.0;
	double vu = 0.0;
	int il = 0;
	int iu = 0;

	double abstol = sqrt(numeric_limits<double>::epsilon()); // absolute tolerance. Shouldn't be too close to epsilon since I don't need very precise estimation of small eigenvalues

	int M   = N;
	int ldz = N;

	vector<int> isuppz(2*M, 0);
	vector<double> work(1, 0.0);        // workspace; size will be determined
	int lwork = -1;          // to start; this lets us determine workspace size
	vector<int> iwork(1, 0); // integer workspace; size to be calculated
	int liwork = -1;         // to start; this lets us determine integer workspace size
	int info = 0;

	dsyevr_(&jobz, &range, &uplo, &N, data_, &lda, &vl, &vu, &il, &iu, &abstol, &M, lam.data(), U.data_, &ldz, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork, &info);

	lwork  = work[0];
	work.resize(static_cast<size_t>(lwork), 0.0);
	liwork = iwork[0];
	iwork.resize(static_cast<size_t>(liwork), 0);

	// run the actual estimation
	dsyevr_(&jobz, &range, &uplo, &N, data_, &lda, &vl, &vu, &il, &iu, &abstol, &M, lam.data(), U.data_, &ldz, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork, &info);

	// set tiny eigenvalues to exactly zero
	for (auto lamIt = lam.begin(); lamIt != lam.end(); ++lamIt) {
		if (fabs(*lamIt) <= abstol) {
			(*lamIt) = 0.0;
		}
	}

}

void Matrix::eigen(const char &tri, const size_t &n, Matrix &U, vector<double> &lam){
#ifndef LMRG_CHECK_OFF
	if (Nrow_ != Ncol_) {
		throw string("ERROR: matrix has to be at least square in eigen()");
	}
	if (Nrow_ < n) {
		throw string("ERROR: the input number of eigenvalues greater than matrix dimensions");
	}
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw string("ERROR: one of the dimensions is zero");
	}
#endif

	if (Nrow_ == n) { // if we are doing all of them, just run regular eigen()
		Matrix::eigen(tri, U, lam);
		return;
	}


	char jobz  = 'V'; // computing eigenvectors
	char range = 'I'; // doing some of them
	char uplo;
	if (tri == 'u') {
		uplo = 'U';
	} else if (tri == 'l'){
		uplo = 'L';
	} else {
		throw string("ERROR: unknown triangle indicator in eigen()");
	}
	int N   = static_cast<int>(Nrow_);
	int lda = static_cast<int>(Nrow_);
	// placeholder variables. Not referenced since we are computing a certain number of eigenvectors, not based on the values of the eigenvalues
	double vl = 0.0;
	double vu = 0.0;
	int il = N - static_cast<int>(n) + 1; // looks like the count base-1
	int iu = N; // do all the remaining eigenvalues

	double abstol = sqrt(numeric_limits<double>::epsilon()); // absolute tolerance. Shouldn't be too close to epsilon since I don't need very precise estimation of small eigenvalues

	int M   = iu - il + 1;
	int ldz = N;

	// test the output size and adjust if necessary
	if ((Nrow_ > U.Nrow_) || (static_cast<size_t>(M) > U.Ncol_)) {
		U.resize(Ncol_, static_cast<size_t>(M));
	}
	if (static_cast<size_t>(M) > lam.size()) {
		lam.resize(static_cast<size_t>(M), 0.0);
	}

	vector<int> isuppz(2*M, 0);
	vector<double> work(1, 0.0);        // workspace; size will be determined
	int lwork = -1;          // to start; this lets us determine workspace size
	vector<int> iwork(1, 0); // integer workspace; size to be calculated
	int liwork = -1;         // to start; this lets us determine integer workspace size
	int info = 0;

	dsyevr_(&jobz, &range, &uplo, &N, data_, &lda, &vl, &vu, &il, &iu, &abstol, &M, lam.data(), U.data_, &ldz, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork, &info);

	lwork  = work[0];
	work.resize(static_cast<size_t>(lwork), 0.0);
	liwork = iwork[0];
	iwork.resize(static_cast<size_t>(liwork), 0);

	// run the actual estimation
	dsyevr_(&jobz, &range, &uplo, &N, data_, &lda, &vl, &vu, &il, &iu, &abstol, &M, lam.data(), U.data_, &ldz, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork, &info);

	// set tiny eigenvalues to exactly zero
	for (auto lamIt = lam.begin(); lamIt != lam.end(); ++lamIt) {
		if (fabs(*lamIt) <= abstol) {
			(*lamIt) = 0.0;
		}
	}

}

void Matrix::eigenSafe(const char &tri, Matrix &U, vector<double> &lam){
#ifndef LMRG_CHECK_OFF
	if (Nrow_ != Ncol_) {
		throw("ERROR: matrix has to be at least square in eigen()");
	}
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw("ERROR: one of the dimensions is zero");
	}
#endif

	// test the output size and adjust if necessary
	if ((Ncol_ > U.Nrow_) || (Ncol_ > U.Ncol_)) {
		U.resize(Ncol_, Ncol_);
	}
	if (Ncol_ > lam.size()) {
		lam.resize(Ncol_, 0.0);
	}

	char jobz  = 'V'; // computing eigenvectors
	char range = 'A'; // doing all of them
	char uplo;
	if (tri == 'u') {
		uplo = 'U';
	} else if (tri == 'l'){
		uplo = 'L';
	} else {
		throw("ERROR: unknown triangle indicator in eigen()");
	}
	int N   = static_cast<int>(Nrow_);
	int lda = static_cast<int>(Nrow_);
	// placeholder variables. Not referenced since we are computing all eigenvectors
	double vl = 0.0;
	double vu = 0.0;
	int il = 0;
	int iu = 0;

	double abstol = sqrt(numeric_limits<double>::epsilon()); // absolute tolerance. Shouldn't be too close to epsilon since I don't need very precise estimation of small eigenvalues

	int M   = N;
	int ldz = N;

	vector<int> isuppz(2*M, 0);
	vector<double> work(1, 0.0);        // workspace; size will be determined
	int lwork = -1;          // to start; this lets us determine workspace size
	vector<int> iwork(1, 0); // integer workspace; size to be calculated
	int liwork = -1;         // to start; this lets us determine integer workspace size
	int info = 0;

	double *dataCopy = new double[Nrow_ * Ncol_];
	memcpy(dataCopy, data_, (Nrow_ * Ncol_)*sizeof(double));

	dsyevr_(&jobz, &range, &uplo, &N, data_, &lda, &vl, &vu, &il, &iu, &abstol, &M, lam.data(), U.data_, &ldz, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork, &info);

	lwork  = work[0];
	work.resize(static_cast<size_t>(lwork), 0.0);
	liwork = iwork[0];
	iwork.resize(static_cast<size_t>(liwork), 0);

	// run the actual estimation
	dsyevr_(&jobz, &range, &uplo, &N, data_, &lda, &vl, &vu, &il, &iu, &abstol, &M, lam.data(), U.data_, &ldz, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork, &info);

	memcpy(data_, dataCopy, (Nrow_ * Ncol_)*sizeof(double));
	delete [] dataCopy;

	// set tiny eigenvalues to exactly zero
	for (auto lamIt = lam.begin(); lamIt != lam.end(); ++lamIt) {
		if (fabs(*lamIt) <= abstol) {
			(*lamIt) = 0.0;
		}
	}

}

void Matrix::premultZ(const Matrix &Z){
#ifndef LMRG_CHECK_OFF
	if (Z.getNcols() != Nrow_) {
		throw("ERROR: Incompatible dimensions between Z and M in premultZ()");
	}
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw("ERROR: one of the dimensions is zero");
		exit(24);
	}
#endif
	// build the vector that stores all the new row IDs for each old row (represented by columns of Z)
	vector< vector<size_t> > fac(Z.getNcols());

	for (size_t oldRow = 0; oldRow < Z.getNcols(); oldRow++) {
		for (size_t newRow = 0; newRow < Z.getNrows(); newRow++) {
			if (Z.getElem(newRow, oldRow) == 1.0) {
				fac[oldRow].push_back(newRow);
			} else if (Z.getElem(newRow, oldRow) != 0.0) {
				throw("ERROR: design matrix can only have elements 1 or 0 in premultZ()");
			}
		}
	}

	double *dataCopy = new double[Nrow_ * Ncol_];
	memcpy(dataCopy, data_, (Nrow_ * Ncol_)*sizeof(double));
	delete [] data_;
	data_ = new double[Z.getNrows() * Ncol_];

	for (size_t oldRow = 0; oldRow < Z.getNcols(); oldRow++) {
		// going through all the rows of Z that correspond to the old row of M
		for (auto facIt = fac[oldRow].begin(); facIt != fac[oldRow].end(); ++facIt) {
			// copying the row of M
			for (size_t jCol = 0; jCol < Ncol_; jCol++) {
				data_[Z.getNrows()*jCol + (*facIt)] = dataCopy[Nrow_*jCol + oldRow];
			}
		}

	}
	Nrow_ = Z.getNrows();

	delete [] dataCopy;
}

void Matrix::premultZ(const Matrix &Z, Matrix &out) const {
#ifndef LMRG_CHECK_OFF
	if (Z.getNcols() != Nrow_) {
		throw("ERROR: Incompatible dimensions between Z and M in premultZ()");
	}
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw("ERROR: one of the dimensions is zero");
	}
#endif
	// build the vector that stores all the new row IDs for each old row (represented by columns of Z)
	vector< vector<size_t> > fac(Z.getNcols());

	for (size_t oldRow = 0; oldRow < Z.getNcols(); oldRow++) {
		for (size_t newRow = 0; newRow < Z.getNrows(); newRow++) {
			if (Z.getElem(newRow, oldRow) == 1.0) {
				fac[oldRow].push_back(newRow);
			} else if (Z.getElem(newRow, oldRow) != 0.0) {
				throw("ERROR: design matrix can only have elements 1 or 0 in premultZ()");
			}
		}
	}

	if ((Z.getNrows() != out.getNrows()) || (Ncol_ != out.getNcols())) {
		out.resize(Z.getNrows(), Ncol_);
	}

	for (size_t oldRow = 0; oldRow < Z.getNcols(); oldRow++) {
		// going through all the rows of Z that correspond to the old row of M
		for (auto facIt = fac[oldRow].begin(); facIt != fac[oldRow].end(); ++facIt) {
			// copying the row of M
			for (size_t jCol = 0; jCol < Ncol_; jCol++) {
				out.data_[Z.getNrows()*jCol + (*facIt)] = data_[Nrow_*jCol + oldRow];
			}
		}

	}

}

void Matrix::premultZt(const Matrix &Z){
#ifndef LMRG_CHECK_OFF
	if (Z.getNrows() != Nrow_) {
		throw("ERROR: Incompatible dimensions between Z and M in premultZt()");
	}
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw("ERROR: one of the dimensions is zero");
	}
#endif
	// build the vector that stores all the new row IDs for each old row (represented by columns of Z)
	vector< vector<size_t> > fac(Z.getNcols());

	for (size_t oldRow = 0; oldRow < Z.getNcols(); oldRow++) {
		for (size_t newRow = 0; newRow < Z.getNrows(); newRow++) {
			if (Z.getElem(newRow, oldRow) == 1.0) {
				fac[oldRow].push_back(newRow);
			} else if (Z.getElem(newRow, oldRow) != 0.0) {
				throw("ERROR: design matrix can only have elements 1 or 0 in premultZt()");
			}
		}
	}

	double *dataCopy = new double[Nrow_ * Ncol_];
	memcpy(dataCopy, data_, (Nrow_ * Ncol_)*sizeof(double));
	delete [] data_;
	data_ = new double[Z.getNcols() * Ncol_](); // value-initializing for summing to work

	for (size_t newRow = 0; newRow < Z.getNcols(); newRow++) {
		// going through all the rows of Z that correspond to the new (summed) row of M
		for (auto facIt = fac[newRow].begin(); facIt != fac[newRow].end(); ++facIt) {
			// summing the rows of M within the group defined by rows of Z
			for (size_t jCol = 0; jCol < Ncol_; jCol++) {
				data_[Z.getNcols()*jCol + newRow] += dataCopy[Nrow_*jCol + (*facIt)];
			}
		}

	}
	Nrow_ = Z.getNcols();

	delete [] dataCopy;

}

void Matrix::premultZt(const Matrix &Z, Matrix &out) const {
#ifndef LMRG_CHECK_OFF
	if (Z.getNrows() != Nrow_) {
		throw("ERROR: Incompatible dimensions between Z and M in premultZt()");
	}
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw("ERROR: one of the dimensions is zero");
	}
#endif
	// build the vector that stores all the new row IDs for each old row (represented by columns of Z)
	vector< vector<size_t> > fac(Z.getNcols());

	for (size_t oldRow = 0; oldRow < Z.getNcols(); oldRow++) {
		for (size_t newRow = 0; newRow < Z.getNrows(); newRow++) {
			if (Z.getElem(newRow, oldRow) == 1.0) {
				fac[oldRow].push_back(newRow);
			} else if (Z.getElem(newRow, oldRow) != 0.0) {
				throw("ERROR: design matrix can only have elements 1 or 0 in premultZt()");
			}
		}
	}

	if ((Z.getNcols() != out.getNrows()) || (Ncol_ != out.getNcols())) {
		out.resize(Z.getNcols(), Ncol_); // resizing already sets all to 0.0
	} else {
		fill(out.data_, out.data_ + (out.Ncol_*out.Nrow_), 0.0);
	}

	for (size_t newRow = 0; newRow < Z.getNcols(); newRow++) {
		// going through all the rows of Z that correspond to the new row of M
		for (auto facIt = fac[newRow].begin(); facIt != fac[newRow].end(); ++facIt) {
			// summing the rows of M within the group defined by rows of Z
			for (size_t jCol = 0; jCol < Ncol_; jCol++) {
				out.data_[Z.getNcols()*jCol + newRow] += data_[Nrow_*jCol + (*facIt)];
			}
		}

	}

}

void Matrix::postmultZ(const Matrix &Z){
#ifndef LMRG_CHECK_OFF
	if (Z.getNrows() != Ncol_) {
		throw("ERROR: Incompatible dimensions between Z and M in postmultZ()");
	}
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw("ERROR: one of the dimensions is zero");
	}
#endif
	// build the vector that stores all the old row IDs for each new row
	vector< vector<size_t> > fac(Z.getNcols());

	for (size_t newCol = 0; newCol < Z.getNcols(); newCol++) {
		for (size_t oldRow = 0; oldRow < Z.getNrows(); oldRow++) {
			if (Z.getElem(oldRow, newCol) == 1.0) {
				fac[newCol].push_back(oldRow);
			} else if (Z.getElem(oldRow, newCol) != 0.0) {
				throw("ERROR: design matrix can only have elements 1 or 0 in postmultZ()");
			}
		}
	}

	double *dataCopy = new double[Nrow_ * Ncol_];
	memcpy(dataCopy, data_, (Nrow_ * Ncol_)*sizeof(double));
	delete [] data_;
	data_ = new double[Nrow_ * Z.getNcols()](); // value initialization for sums to work

	for (size_t newCol = 0; newCol < Z.getNcols(); newCol++) {
		// going through all the rows of Z that correspond to the new column of M
		for (auto facIt = fac[newCol].begin(); facIt != fac[newCol].end(); ++facIt) {
			// summing the rows of M within the group defined by rows of Z
			for (size_t iRow = 0; iRow < Nrow_; iRow++) {
				data_[Z.getNcols()*newCol + iRow] += dataCopy[Nrow_*(*facIt) + iRow];
			}
		}

	}
	Ncol_ = Z.getNcols();
	delete [] dataCopy;
}


void Matrix::postmultZ(const Matrix &Z, Matrix &out) const{
#ifndef LMRG_CHECK_OFF
	if (Z.getNrows() != Ncol_) {
		throw("ERROR: Incompatible dimensions between Z and M in postmultZ()");
	}
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw("ERROR: one of the dimensions is zero");
	}
#endif
	// build the vector that stores all the old column IDs for each new column
	vector< vector<size_t> > fac(Z.getNcols());

	for (size_t newCol = 0; newCol < Z.getNcols(); newCol++) {
		for (size_t oldRow = 0; oldRow < Z.getNrows(); oldRow++) {
			if (Z.getElem(oldRow, newCol) == 1.0) {
				fac[newCol].push_back(oldRow);
			} else if (Z.getElem(oldRow, newCol) != 0.0) {
				throw("ERROR: design matrix can only have elements 1 or 0 in postmultZ()");
			}
		}
	}

	if ((Z.getNcols() != out.getNcols()) || (Nrow_ != out.getNrows())) {
		out.resize(Nrow_, Z.getNcols());
	}

	for (size_t newCol = 0; newCol < Z.getNcols(); newCol++) {
		// going through all the rows of Z that correspond to the new column of M
		for (auto facIt = fac[newCol].begin(); facIt != fac[newCol].end(); ++facIt) {
			// summing the rows of M within the group defined by rows of Z
			for (size_t iRow = 0; iRow < Nrow_; iRow++) {
				out.data_[Z.getNcols()*newCol + iRow] += data_[Nrow_*(*facIt) + iRow];
			}
		}

	}

}
void Matrix::postmultZt(const Matrix &Z){
#ifndef LMRG_CHECK_OFF
	if (Z.getNcols() != Ncol_) { // Z not transposed
		throw("ERROR: Incompatible dimensions between Z and M in postmultZ()");
	}
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw("ERROR: one of the dimensions is zero");
	}
#endif
	// build the vector that stores all the new row IDs for each old row (represented by columns of Z)
	vector< vector<size_t> > fac(Z.getNcols());

	for (size_t oldRow = 0; oldRow < Z.getNcols(); oldRow++) {
		for (size_t newRow = 0; newRow < Z.getNrows(); newRow++) {
			if (Z.getElem(newRow, oldRow) == 1.0) {
				fac[oldRow].push_back(newRow);
			} else if (Z.getElem(newRow, oldRow) != 0.0) {
				throw("ERROR: design matrix can only have elements 1 or 0 in postmultZ()");
			}
		}
	}

	double *dataCopy = new double[Nrow_ * Ncol_];
	memcpy(dataCopy, data_, (Nrow_ * Ncol_)*sizeof(double));
	delete [] data_;
	data_ = new double[Z.getNrows() * Nrow_];

	for (size_t oldCol = 0; oldCol < Z.getNcols(); oldCol++) {
		// going through all the rows of Z that correspond to the old column of M
		for (auto facIt = fac[oldCol].begin(); facIt != fac[oldCol].end(); ++facIt) {
			// copying the column of M
			memcpy(data_ + (*facIt)*Nrow_, dataCopy + oldCol*Nrow_, Nrow_*sizeof(double));
		}
	}

	Ncol_ = Z.getNrows();
	delete [] dataCopy;
}


void Matrix::postmultZt(const Matrix &Z, Matrix &out) const{
#ifndef LMRG_CHECK_OFF
	if (Z.getNcols() != Ncol_) { // Z not transposed
		throw("ERROR: Incompatible dimensions between Z and M in postmultZ()");
	}
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw("ERROR: one of the dimensions is zero");
	}
#endif
	// build the vector that stores all the new row IDs for each old row (represented by columns of Z)
	vector< vector<size_t> > fac(Z.getNcols());

	for (size_t oldRow = 0; oldRow < Z.getNcols(); oldRow++) {
		for (size_t newRow = 0; newRow < Z.getNrows(); newRow++) {
			if (Z.getElem(newRow, oldRow) == 1.0) {
				fac[oldRow].push_back(newRow);
			} else if (Z.getElem(newRow, oldRow) != 0.0) {
				throw("ERROR: design matrix can only have elements 1 or 0 in postmultZ()");
			}
		}
	}

	// careful: Z is not transposed, but will be constructing MZ^t
	if ((Z.getNrows() != out.getNcols()) || (Nrow_ != out.getNrows())) {
		out.resize(Z.getNrows(), Nrow_);
	}

	for (size_t oldCol = 0; oldCol < Z.getNcols(); oldCol++) {
		// going through all the rows of Z that correspond to the old column of M
		for (auto facIt = fac[oldCol].begin(); facIt != fac[oldCol].end(); ++facIt) {
			// copying the column of M
			memcpy(out.data_ + (*facIt)*Nrow_, data_ + oldCol*Nrow_, Nrow_*sizeof(double));
		}
	}

}

void Matrix::syrk(const char &tri, const double &alpha, const double &beta, Matrix &C) const {
#ifndef LMRG_CHECK_OFF
	if ((Ncol_ > INT_MAX) || (Nrow_ > INT_MAX)) {
		throw("ERROR: at least one matrix dimension too big to safely convert to int in syrk()");
	}
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw("ERROR: one of the dimensions is zero");
	}
#endif

	if ((C.getNrows() != Ncol_) || (C.getNcols() != Ncol_)) {
		C.resize(Ncol_, Ncol_);
	}

	// integer parameters
	const int n   = static_cast<int>(Ncol_);
	const int k   = static_cast<int>(Nrow_);
	const int lda = static_cast<int>(Nrow_);
	const int ldc = static_cast<int>(Ncol_);

	// transpose token
	const char trans = 't';

	dsyrk_(&tri, &trans, &n, &k, &alpha, data_, &lda, &beta, C.data_, &ldc);

}

void Matrix::tsyrk(const char &tri, const double &alpha, const double &beta, Matrix &C) const {
#ifndef LMRG_CHECK_OFF
	if ((Ncol_ > INT_MAX) || (Nrow_ > INT_MAX)) {
		throw("ERROR: at least one matrix dimension too big to safely convert to int in tsyrk()");
	}
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw("ERROR: one of the dimensions is zero");
	}
#endif

	if ((C.getNrows() != Nrow_) || (C.getNcols() != Nrow_)) {
		C.resize(Nrow_, Nrow_);
	}

	// integer parameters
	const int n   = static_cast<int>(Nrow_);
	const int k   = static_cast<int>(Ncol_);
	const int lda = static_cast<int>(Nrow_);
	const int ldc = static_cast<int>(Nrow_);

	// transpose token
	const char trans = 'n';

	dsyrk_(&tri, &trans, &n, &k, &alpha, data_, &lda, &beta, C.data_, &ldc);
}

void Matrix::symm(const char &tri, const char &side, const double &alpha, const Matrix &symA, const double &beta, Matrix &C) const{
#ifndef LMRG_CHECK_OFF
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw("ERROR: one of the dimensions is zero");
	}
	if (symA.getNrows() != symA.getNcols()) {
		throw("ERROR: symmetric matrix symA has to be square in symm()");
	}
	if (side == 'l') {
		if ((Nrow_ > INT_MAX) || (symA.getNcols() > INT_MAX)) {
		throw("ERROR: at least one matrix dimension too big to safely convert to int in symm()");
		}
	} else if (side == 'r') {
		if ((symA.getNrows() > INT_MAX) || (Ncol_ > INT_MAX)) {
			throw("ERROR: at least one matrix dimension too big to safely convert to int in symm()");
		}
	}

	if ((symA.getNcols() != Nrow_) && (side == 'l')) { // AB
		throw("ERROR: Incompatible dimensions between B and A in symm()");
	}
	if ((symA.getNrows() != Ncol_) && (side == 'r')) { // BA
		throw("ERROR: Incompatible dimensions between A and B in symm()");
	}
#endif

	int m;
	int n;
	if (side == 'l') { // AB
		m = static_cast<int>(symA.getNrows());
		n = static_cast<int>(Ncol_);
		if ((C.getNrows() != symA.getNrows()) || (C.getNcols() != Ncol_)) {
			C.resize(symA.getNrows(), Ncol_);
		}
	} else if (side == 'r') { // BA
		m = static_cast<int>(Nrow_);
		n = static_cast<int>(symA.getNcols());
		if ((C.getNrows() != Nrow_) || (C.getNcols() != symA.getNcols())) {
			C.resize(Nrow_, symA.getNcols());
		}
	} else {
		throw("ERROR: unknown side indicator in symm()");
	}

	// final integer parameters
	const int lda = static_cast<int>(symA.getNrows());
	const int ldb = static_cast<int>(Nrow_);
	const int ldc = m; // for clarity

	dsymm_(&side, &tri, &m, &n, &alpha, symA.data_, &lda, data_, &ldb, &beta, C.data_, &ldc);

}

void Matrix::symc(const char &tri, const double &alpha, const Matrix &X, const size_t &xCol, const double &beta, vector<double> &y) const{
#ifndef LMRG_CHECK_OFF
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw("ERROR: one of the dimensions is zero");
	}
	if (Ncol_ != Nrow_) {
		throw("ERROR: symmetric matrix (current object) has to be square in symc()");
	}
	if ((Ncol_ > INT_MAX) || (X.getNrows() > INT_MAX)) {
		throw("ERROR: at least one matrix dimension too big to safely convert to int in symc()");
	}
	if (X.getNrows() != Ncol_) {
		throw("ERROR: Incompatible dimensions between A and X in symc()");
	}
	if (xCol >= X.getNcols()) {
		throw("ERROR: column index out of range for matrix X in symc()");
	}
#endif
	if (y.size() < Nrow_) {
		y.resize(Nrow_);
	}

	// BLAS routine constants
	const int n    = static_cast<int>(Nrow_);
	const int lda  = n;
	const int incx = 1;
	const int incy = 1;

	const double *xbeg = X.data_ + xCol*(X.Nrow_); // offset to the column of interest

	dsymv_(&tri, &n, &alpha, data_, &lda, xbeg, &incx, &beta, y.data(), &incy);
}

void Matrix::gemm(const bool &transA, const double &alpha, const Matrix &A, const bool &transB, const double &beta, Matrix &C) const{
#ifndef LMRG_CHECK_OFF
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw("ERROR: one of the dimensions is zero");
	}
	if ((A.getNcols() > INT_MAX) || (A.getNrows() > INT_MAX)) {
		throw("ERROR: at least one A matrix dimension too big to safely convert to int in gemm()");
	}

	if (transB) {
		if (Nrow_ > INT_MAX) {
			throw("ERROR: at least one B matrix dimension too big to safely convert to int in gemm()");
		}
	} else {
		if (Ncol_ > INT_MAX) {
			throw("ERROR: at least one B matrix dimension too big to safely convert to int in gemm()");
		}
	}
	if (transA) {
		if (transB && (A.getNrows() != Ncol_)) {
			throw("ERROR: Incompatible dimensions between A^T and B^T in gemm()");
		} else if (!transB && (A.getNrows() != Nrow_)){
			throw("ERROR: Incompatible dimensions between A^T and B in gemm()");
		}

	} else {
		if (transB && (A.getNcols() != Ncol_)) {
			throw("ERROR: Incompatible dimensions between A and B^T in gemm()");
			exit(16);
		} else if (!transB && (A.getNcols() != Nrow_)) {
			throw("ERROR: Incompatible dimensions between A and B in gemm()");
		}
	}
#endif

	char tAtok;
	char tBtok;

	int m;
	int k;
	int n;
	if (transA) {
		tAtok = 't';
		m     = static_cast<int>(A.getNcols());
		k     = static_cast<int>(A.getNrows());
		if (transB) {
			tBtok = 't';
			n     = static_cast<int>(Nrow_);
			if ((C.getNrows() != A.getNcols()) || (C.getNcols() != Nrow_)) {
				C.resize(A.getNcols(), Nrow_);
			}
		} else {
			tBtok = 'n';
			n     = static_cast<int>(Ncol_);
			if ((C.getNrows() != A.getNcols()) || (C.getNcols() != Ncol_)) {
				C.resize(A.getNcols(), Ncol_);
			}
		}
	} else {
		tAtok = 'n';
		m     = static_cast<int>(A.getNrows());
		k     = static_cast<int>(A.getNcols());
		if (transB) {
			tBtok = 't';
			n     = static_cast<int>(Nrow_);
			if ((C.getNrows() != A.getNrows()) || (C.getNcols() != Nrow_)) {
				C.resize(A.getNrows(), Nrow_);
			}
		} else {
			tBtok = 'n';
			n     = static_cast<int>(Ncol_);
			if ((C.getNrows() != A.getNrows()) || (C.getNcols() != Ncol_)) {
				C.resize(A.getNrows(), Ncol_);
			}
		}
	}

	const int lda = (transA ? k : m);
	const int ldb = (transB ? n : k);
	const int ldc = m;

	dgemm_(&tAtok, &tBtok, &m, &n, &k, &alpha, A.data_, &lda, data_, &ldb, &beta, C.data_, &ldc);

}
void Matrix::gemc(const bool &trans, const double &alpha, const Matrix &X, const size_t &xCol, const double &beta, vector<double> &y) const {
#ifndef LMRG_CHECK_OFF
	if ( (Nrow_ == 0) || (Ncol_ == 0) ) {
		throw("ERROR: one of the dimensions is zero");
	}
	if (trans) {
		if ((Nrow_ > INT_MAX) || (X.getNrows() > INT_MAX)) {
			throw("ERROR: at least one matrix dimension too big to safely convert to int in gemc()");
		}
		if (Nrow_ != X.getNrows()) {
			throw("ERROR: Incompatible dimensions between A and X in gemc()");

		}
	} else {
		if ((Ncol_ > INT_MAX) || (X.getNrows() > INT_MAX)) {
			throw("ERROR: at least one matrix dimension too big to safely convert to int in gemc()");
		}
		if (Ncol_ != X.getNrows()) {
			throw("ERROR: Incompatible dimensions between A and X in gemc()");

		}
	}
	if (xCol >= X.getNcols()) {
		throw("ERROR: column index out of range for matrix X in gemc()");
	}
#endif

	if (y.size() < Nrow_) {
		y.resize(Nrow_);
	}

	// Establish constants for DGEMV
	const char tTok = (trans ? 't' : 'n');

	const int m    = static_cast<int>(Nrow_);
	const int n    = static_cast<int>(Ncol_);
	const int lda  = m;
	const int incx = 1;
	const int incy = 1;

	const double *xbeg = X.data_ + xCol*(X.Nrow_); // offset to the column of interest

	dgemv_(&tTok, &m, &n, &alpha, data_, &lda, xbeg, &incx, &beta, y.data(), &incy);

}

Matrix Matrix::colShuffle() const{
	Matrix out(Nrow_, Ncol_);
	RanDraw rng;
	vector<uint64_t> ind = rng.shuffleUint(Ncol_); // shuffled column index
	for (size_t j = 0; j < Ncol_; ++j) {
		for (size_t i = 0; i < Nrow_; ++i) {
			out.setElem(i, ind[j], this->getElem(i, j));
		}
	}

	return out;
}

Matrix Matrix::rowShuffle() const{
	Matrix out(Nrow_, Ncol_);
	RanDraw rng;
	vector<uint64_t> ind = rng.shuffleUint(Nrow_); // shuffled column index
	for (size_t j = 0; j < Ncol_; ++j) {
		for (size_t i = 0; i < Nrow_; ++i) {
			out.setElem(ind[i], j, this->getElem(i, j));
		}
	}

	return out;
}

Matrix Matrix::operator*(const Matrix &m) const{
#ifndef LMRG_CHECK_OFF
	if ((Nrow_ != m.Nrow_) || (Ncol_ != m.Ncol_)) {
		throw("ERROR: Incompatible dimensions between matrices in the Hadamard product");
	}
#endif
	Matrix res(*this);
	for (size_t iElm = 0; iElm < Ncol_*Nrow_; iElm++) {
		res.data_[iElm] *= m.data_[iElm];
	}

	return res;
}

Matrix Matrix::operator*(const double &scal) const{
	Matrix res(*this);
	for (size_t iElm = 0; iElm < Ncol_*Nrow_; iElm++) {
		res.data_[iElm] *= scal;
	}

	return res;
}

Matrix Matrix::operator/(const Matrix &m) const{
#ifndef LMRG_CHECK_OFF
	if ((Nrow_ != m.Nrow_) || (Ncol_ != m.Ncol_)) {
		throw("ERROR: Incompatible dimensions between matrices in the Hadamard product");
	}
#endif
	Matrix res(*this);
	for (size_t iElm = 0; iElm < Ncol_*Nrow_; iElm++) {
		res.data_[iElm] /= m.data_[iElm];
	}

	return res;
}

Matrix Matrix::operator/(const double &scal) const{
	Matrix res(*this);
	for (size_t iElm = 0; iElm < Ncol_*Nrow_; iElm++) {
		res.data_[iElm] /= scal;
	}

	return res;
}

Matrix Matrix::operator+(const Matrix &m) const{
#ifndef LMRG_CHECK_OFF
	if ((Nrow_ != m.Nrow_) || (Ncol_ != m.Ncol_)) {
		throw("ERROR: Incompatible dimensions between matrices in the Hadamard product");
	}
#endif
	Matrix res(*this);
	for (size_t iElm = 0; iElm < Ncol_*Nrow_; iElm++) {
		res.data_[iElm] += m.data_[iElm];
	}

	return res;
}

Matrix Matrix::operator+(const double &scal) const{
	Matrix res(*this);
	for (size_t iElm = 0; iElm < Ncol_*Nrow_; iElm++) {
		res.data_[iElm] += scal;
	}

	return res;
}
Matrix Matrix::operator-(const Matrix &m) const{
#ifndef LMRG_CHECK_OFF
	if ((Nrow_ != m.Nrow_) || (Ncol_ != m.Ncol_)) {
		throw("ERROR: Incompatible dimensions between matrices in the Hadamard product");
	}
#endif
	Matrix res(*this);
	for (size_t iElm = 0; iElm < Ncol_*Nrow_; iElm++) {
		res.data_[iElm] -= m.data_[iElm];
	}

	return res;
}

Matrix Matrix::operator-(const double &scal) const{
	Matrix res(*this);
	for (size_t iElm = 0; iElm < Ncol_*Nrow_; iElm++) {
		res.data_[iElm] -= scal;
	}

	return res;
}

Matrix& Matrix::operator+=(const double &scal){
	for (size_t iElm = 0; iElm < Ncol_*Nrow_; iElm++) {
		data_[iElm] += scal;
	}

	return *this;
}

Matrix& Matrix::operator*=(const double &scal){
	for (size_t iElm = 0; iElm < Ncol_*Nrow_; iElm++) {
		data_[iElm] *= scal;
	}

	return *this;
}

Matrix& Matrix::operator-=(const double &scal){
	for (size_t iElm = 0; iElm < Ncol_*Nrow_; iElm++) {
		data_[iElm] -= scal;
	}

	return *this;
}

Matrix& Matrix::operator/=(const double &scal){
	for (size_t iElm = 0; iElm < Ncol_*Nrow_; iElm++) {
		data_[iElm] /= scal;
	}

	return *this;
}

void Matrix::rowMeans(vector<double> &means) const{
	if (means.size() < Nrow_) {
		means.resize(Nrow_);
	}
	for (size_t iRow = 0; iRow < Nrow_; iRow++) {
		means[iRow] = 0.0; // in case something was in the vector passed to the function and resize did not erase it
		for (size_t jCol = 0; jCol < Ncol_; jCol++) {
			// numerically stable recursive mean calculation. GSL does it this way.
			means[iRow] += (data_[Nrow_*jCol + iRow] - means[iRow])/static_cast<double>(jCol + 1);
		}
	}
}
void Matrix::colMeans(vector<double> &means) const{
	if (means.size() < Ncol_) {
		means.resize(Ncol_);
	}
	for (size_t jCol = 0; jCol < Ncol_; jCol++) {
		means[jCol] = 0.0; // in case something was in the vector passed to the function and resize did not erase it
		for (size_t iRow = 0; iRow < Nrow_; iRow++) {
			// numerically stable recursive mean calculation. GSL does it this way.
			means[jCol] += (data_[Nrow_*jCol + iRow] - means[jCol])/static_cast<double>(iRow + 1);
		}
	}
}
void Matrix::rowSums(vector<double> &sums) const{
	if (sums.size() < Nrow_) {
		sums.resize(Nrow_);
	}
	for (size_t iRow = 0; iRow < Nrow_; iRow++) {
		sums[iRow] = 0.0; // in case something was in the vector passed to the function and resize did not erase it
		for (size_t jCol = 0; jCol < Ncol_; jCol++) {
			// not necessarily mumerically stable. Revisit later
			sums[iRow] += data_[Nrow_*jCol + iRow];
		}
	}
}
void Matrix::colSums(vector<double> &sums) const{
	if (sums.size() < Ncol_) {
		sums.resize(Ncol_);
	}
	for (size_t jCol = 0; jCol < Ncol_; jCol++) {
		sums[jCol] = 0.0; // in case something was in the vector passed to the function and resize did not erase it
		for (size_t iRow = 0; iRow < Nrow_; iRow++) {
			// not necessarily mumerically stable. Revisit later
			sums[jCol] += data_[Nrow_*jCol + iRow];
		}
	}
}
void Matrix::rowMultiply(const vector<double> &scalars){
#ifndef LMRG_CHECK_OFF
	if (scalars.size() != Ncol_) {
		throw("ERROR: Vector of scalars has wrong length in rowMultiply(vector)");
	}
#endif

	for (size_t iRow = 0; iRow < Nrow_; iRow++) {
		for (size_t jCol = 0; jCol < Ncol_; jCol++) {
			data_[Nrow_*jCol + iRow] *= scalars[jCol];
		}
	}
}
void Matrix::rowMultiply(const double &scalar, const size_t &iRow){
#ifndef LMRG_CHECK_OFF
	if (iRow >= Nrow_) {
		throw("ERROR: Row index out of bounds in rowMultiply(scalar)");
	}
#endif
	for (size_t jCol = 0; jCol < Ncol_; jCol++) {
		data_[Nrow_*jCol + iRow] *= scalar;
	}
}
void Matrix::colMultiply(const vector<double> &scalars){
#ifndef LMRG_CHECK_OFF
	if (scalars.size() != Nrow_) {
		throw("ERROR: Vector of scalars has wrong length in colMultiply(vector)");
	}
#endif

	for (size_t jCol = 0; jCol < Ncol_; jCol++) {
		for (size_t iRow = 0; iRow < Nrow_; iRow++) {
			data_[Nrow_*jCol + iRow] *= scalars[iRow];
		}
	}
}
void Matrix::colMultiply(const double &scalar, const size_t &jCol){
#ifndef LMRG_CHECK_OFF
	if (jCol >= Ncol_) {
		throw("ERROR: Column index out of bounds in colMultiply(scalar)");
	}
#endif
	for (size_t iRow = 0; iRow < Nrow_; iRow++) {
		data_[Nrow_*jCol + iRow] *= scalar;
	}
}
void Matrix::rowDivide(const vector<double> &scalars){
#ifndef LMRG_CHECK_OFF
	if (scalars.size() != Ncol_) {
		throw("ERROR: Vector of scalars has wrong length in rowDivide(vector)");
	}
#endif

	for (size_t iRow = 0; iRow < Nrow_; iRow++) {
		for (size_t jCol = 0; jCol < Ncol_; jCol++) {
			data_[Nrow_*jCol + iRow] /= scalars[jCol];
		}
	}
}
void Matrix::rowDivide(const double &scalar, const size_t &iRow){
#ifndef LMRG_CHECK_OFF
	if (iRow >= Nrow_) {
		throw("ERROR: Row index out of bounds in rowDivide(scalar)");
	}
#endif
	for (size_t jCol = 0; jCol < Ncol_; jCol++) {
		data_[Nrow_*jCol + iRow] /= scalar;
	}
}
void Matrix::colDivide(const vector<double> &scalars){
#ifndef LMRG_CHECK_OFF
	if (scalars.size() != Nrow_) {
		throw("ERROR: Vector of scalars has wrong length in colDivide(vector)");
	}
#endif

	for (size_t jCol = 0; jCol < Ncol_; jCol++) {
		for (size_t iRow = 0; iRow < Nrow_; iRow++) {
			data_[Nrow_*jCol + iRow] /= scalars[iRow];
		}
	}
}
void Matrix::colDivide(const double &scalar, const size_t &jCol){
#ifndef LMRG_CHECK_OFF
	if (jCol >= Ncol_) {
		throw("ERROR: Column index out of bounds in colDivide(scalar)");
	}
#endif
	for (size_t iRow = 0; iRow < Nrow_; iRow++) {
		data_[Nrow_*jCol + iRow] /= scalar;
	}
}
void Matrix::rowAdd(const vector<double> &scalars){
#ifndef LMRG_CHECK_OFF
	if (scalars.size() != Ncol_) {
		throw("ERROR: Vector of scalars has wrong length in rowAdd(vector)");
	}
#endif

	for (size_t iRow = 0; iRow < Nrow_; iRow++) {
		for (size_t jCol = 0; jCol < Ncol_; jCol++) {
			data_[Nrow_*jCol + iRow] += scalars[jCol];
		}
	}
}
void Matrix::rowAdd(const double &scalar, const size_t &iRow){
#ifndef LMRG_CHECK_OFF
	if (iRow >= Nrow_) {
		throw("ERROR: Row index out of bounds in rowAdd(scalar)");
	}
#endif
	for (size_t jCol = 0; jCol < Ncol_; jCol++) {
		data_[Nrow_*jCol + iRow] += scalar;
	}
}
void Matrix::colAdd(const vector<double> &scalars){
#ifndef LMRG_CHECK_OFF
	if (scalars.size() != Nrow_) {
		throw("ERROR: Vector of scalars has wrong length in colAdd(vector)");
	}
#endif

	for (size_t jCol = 0; jCol < Ncol_; jCol++) {
		for (size_t iRow = 0; iRow < Nrow_; iRow++) {
			data_[Nrow_*jCol + iRow] += scalars[iRow];
		}
	}
}
void Matrix::colAdd(const double &scalar, const size_t &jCol){
#ifndef LMRG_CHECK_OFF
	if (jCol >= Ncol_) {
		throw("ERROR: Column index out of bounds in colAdd(scalar)");
	}
#endif
	for (size_t iRow = 0; iRow < Nrow_; iRow++) {
		data_[Nrow_*jCol + iRow] += scalar;
	}
}
void Matrix::rowSub(const vector<double> &scalars){
#ifndef LMRG_CHECK_OFF
	if (scalars.size() != Ncol_) {
		throw("ERROR: Vector of scalars has wrong length in rowSub(vector)");
	}
#endif

	for (size_t iRow = 0; iRow < Nrow_; iRow++) {
		for (size_t jCol = 0; jCol < Ncol_; jCol++) {
			data_[Nrow_*jCol + iRow] -= scalars[jCol];
		}
	}
}
void Matrix::rowSub(const double &scalar, const size_t &iRow){
#ifndef LMRG_CHECK_OFF
	if (iRow >= Nrow_) {
		throw("ERROR: Row index out of bounds in rowSub(scalar)");
	}
#endif
	for (size_t jCol = 0; jCol < Ncol_; jCol++) {
		data_[Nrow_*jCol + iRow] -= scalar;
	}
}
void Matrix::colSub(const vector<double> &scalars){
#ifndef LMRG_CHECK_OFF
	if (scalars.size() != Nrow_) {
		throw("ERROR: Vector of scalars has wrong length in colSub(vector)");
	}
#endif

	for (size_t jCol = 0; jCol < Ncol_; jCol++) {
		for (size_t iRow = 0; iRow < Nrow_; iRow++) {
			data_[Nrow_*jCol + iRow] -= scalars[iRow];
		}
	}
}
void Matrix::colSub(const double &scalar, const size_t &jCol){
#ifndef LMRG_CHECK_OFF
	if (jCol >= Ncol_) {
		throw("ERROR: Column index out of bounds in colSub(scalar)");
	}
#endif
	for (size_t iRow = 0; iRow < Nrow_; iRow++) {
		data_[Nrow_*jCol + iRow] -= scalar;
	}
}

void Matrix::appendCol(const Matrix &cols){
#ifndef LMRG_CHECK_OFF
	if (Ncol_ > cols.Ncol_ + Ncol_) { // happens on wrap-around
		throw("ERROR: Number of columns too big to expand");
	}
#endif
	// if the other matrix is empty, do nothing
	if ( (cols.Ncol_ == 0) || (cols.Ncol_ == 0) ) {
		return;
	}
	if (this == &cols) { // self-appending
		if (Ncol_ && Nrow_) {
			double *dataCopy = new double[Nrow_ * Ncol_];
			memcpy(dataCopy, data_, (Nrow_ * Ncol_)*sizeof(double));

			delete [] data_;
			data_ = new double[Nrow_ * 2 * Ncol_];
			memcpy(data_, dataCopy, (Nrow_ * Ncol_)*sizeof(double));
			memcpy(data_ + Nrow_*Ncol_, dataCopy, (Nrow_ * Ncol_)*sizeof(double));
			Ncol_ *= 2;
			delete [] dataCopy;
		} else {
			return;
		}
	} else {
#ifndef LMRG_CHECK_OFF
		if (Nrow_ != cols.Nrow_) {
			throw("ERROR: Number of rows in appeneded object not equal to the number of rows in focal matrix");
		}
#endif
		double *dataCopy = new double[Nrow_ * Ncol_];
		memcpy(dataCopy, data_, (Nrow_ * Ncol_)*sizeof(double));
		delete [] data_;
		data_ = new double[Nrow_ * (cols.Ncol_ + Ncol_)];
		memcpy(data_, dataCopy, (Nrow_ * Ncol_)*sizeof(double));
		delete [] dataCopy;

		memcpy(data_ + Nrow_*Ncol_, cols.data_, (Nrow_ * (cols.Ncol_))*sizeof(double));
		Ncol_ += cols.Ncol_;


	}
}

void Matrix::appendRow(const Matrix &rows){
#ifndef LMRG_CHECK_OFF
	if (Nrow_ > rows.Nrow_ + Nrow_) { // happens on wrap-around
		throw("ERROR: Number of rows too big to expand");
	}
#endif
	// if the other matrix is empty, do nothing
	if ( (rows.Ncol_ == 0) || (rows.Ncol_ == 0) ) {
		return;
	}
	if (this == &rows) { // self-appending
		if (Nrow_ && Ncol_) {
			double *dataCopy = new double[Nrow_ * Ncol_];
			memcpy(dataCopy, data_, (Nrow_ * Ncol_)*sizeof(double));
			delete [] data_;

			data_ = new double[2 * Nrow_ * Ncol_];
			// since rows are discontinuous, have to go column by column and copy rows within each column
			for (size_t jCol = 0; jCol < Ncol_; jCol++) {
				memcpy(data_ + 2*Nrow_*jCol, dataCopy + Nrow_*jCol, Nrow_*sizeof(double));
				memcpy(data_ + 2*Nrow_*jCol + Nrow_, dataCopy + Nrow_*jCol, Nrow_*sizeof(double));
			}
			Nrow_ *= 2;
			delete [] dataCopy;
		} else {
			return;
		}
	} else {
#ifndef LMRG_CHECK_OFF
		if (Ncol_ != rows.Ncol_) {
			throw("ERROR: Number of columns in appeneded object not equal to the number of columns in focal matrix");
		}
#endif
		double *dataCopy = new double[Nrow_ * Ncol_];
		memcpy(dataCopy, data_, (Nrow_ * Ncol_)*sizeof(double));
		delete [] data_;
		data_ = new double[(rows.Nrow_ + Nrow_) * Ncol_];

		// since rows are discontinuous, have to go column by column and copy rows within each column
		for (size_t jCol = 0; jCol < Ncol_; jCol++) {
			memcpy(data_ + (rows.Nrow_+Nrow_)*jCol, dataCopy + Nrow_*jCol, Nrow_*sizeof(double));
			memcpy(data_ + (rows.Nrow_+Nrow_)*jCol + Nrow_, rows.data_ + (rows.Nrow_)*jCol, (rows.Nrow_)*sizeof(double));
		}
		Nrow_ += rows.Nrow_;
		delete [] dataCopy;
	}
}


