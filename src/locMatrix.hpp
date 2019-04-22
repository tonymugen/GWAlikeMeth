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
 * This is the project header file containing class definitions and interface documentation.
 *
 */

#ifndef locMatrix_hpp
#define locMatrix_hpp

#include <vector>
#include <string>

using std::vector;
using std::string;

namespace BayesicSpace {
	class Matrix;
	/** \defgroup arithmetic Arithmetic operators
	 *
	 * Scalar by matrix addition and multiplication operators to maintain commutativity.
	 *
	 * @{
	 */
	/** \brief Scalar-matrix product
	 *
	 * \param[in] scal scalar
	 * \param[in] m matrix
	 * \return Matrix result
	 *
	 */
	Matrix operator*(const double &scal, const Matrix &m);
	/** \brief Scalar-matrix addition
	 *
	 * \param[in] scal scalar
	 * \param[in] m matrix
	 * \return Matrix result
	 *
	 */
	Matrix operator+(const double &scal, const Matrix &m);
	/** @} */

	/** \brief Test matrix class
	 *
	 * A test matrix class. The data type is _double_. The matrix is column-major to comply with LAPACK and BLAS routines.
	 * Columns and rows are base-0. Range checking is done unless the flag -DLMRG_CHECK_OFF is set at compile time.
	 *
	 * TODO: transpose
	 * TODO: make thread-safe with C++-11 facilities
	 *
	 */
	class Matrix {
		friend Matrix operator+(const double &scal, const Matrix &m);
		friend Matrix operator*(const double &scal, const Matrix &m);

	public:
		/** \brief Default constructor
		 *
		 */
		Matrix() : data_{nullptr}, Nrow_ (0), Ncol_ (0) {}; // setting to nullptr is important for safety of delete [] when compiled with -O3
		/** \brief Allocating constructor
		 *
		 * Allocates memory but does not initialize
		 *
		 * \param[in] nrow number of rows
		 * \param[in] ncol number of columns
		 *
		 */
		Matrix(const size_t &nrow, const size_t &ncol);
		/** \brief Initializing constructor
		 *
		 * Allocates memory and initializes the elements to the given value.
		 *
		 * \param[in] val initializing value
		 * \param[in] nrow number of rows
		 * \param[in] ncol number of columns
		 *
		 */
		Matrix(const double &val, const size_t &nrow, const size_t &ncol);
		/** \brief Constructor from C array
		 *
		 * Allocates memory and initializes by copying values from the given array. The user is resposible for making sure the input array is long enough.
		 *
		 * \param[in] inArr initializing array
		 * \param[in] nrow number of rows
		 * \param[in] ncol number of columns
		 *
		 */
		Matrix(const double inArr[], const size_t &nrow, const size_t &ncol);
		/** \brief Constructor from C++ vector
		 *
		 * Allocates memory and initializes by copying values from the given vector. Vector can be bigger (but not smaller) than necessary. If so, the beginning \f$N_{row} \times N_{col}\f$ elements are used.
		 *
		 * \param[in] inVec initializing array
		 * \param[in] nrow number of rows
		 * \param[in] ncol number of columns
		 */
		Matrix(const vector<double> &inVec, const size_t &nrow, const size_t &ncol);
		/** \brief Constructor from file
		 *
		 *	Gets data from a space-delimited file.
		 *
		 * \param[in] fileName input file name
		 * \param[in] delim column delimiter
		 *
		 */
		Matrix(const string &fileName, const char &delim);

		/** \brief Destructor
		 *
		 */
		~Matrix();

		/** \brief Copy constructor
		 *
		 * \param[in] inMat object to be copied
		 */
		Matrix(const Matrix &inMat);

		/** \brief Copy assignment operator
		 *
		 * \param[in] inMat object to be copied
		 * \return Matrix target object
		 *
		 */
		Matrix& operator=(const Matrix &inMat);
		/** \brief Move constructor
		 *
		 * \param[in] inMat object to be moved
		 */
		Matrix(Matrix &&inMat);
		/** \brief Move assignment operator
		 *
		 * \param[in] inMat object to be moved
		 * \return Matrix target object
		 *
		 */
		Matrix& operator=(Matrix &&inMat);

		/** \brief Access to number of rows
		 *
		 * \return size_t number of rows
		 */
		size_t getNrows() const{return Nrow_; };
		/** \brief Access to number of columns
		 *
		 * \return size_t number of columns
		 */
		size_t getNcols() const{return Ncol_; };
		/** \brief Access to an element
		 *
		 * \param[in] iRow row number
		 * \param[in] jCol column number
		 * \return double element value
		 */
		double getElem(const size_t &iRow, const size_t &jCol) const;
		/** \brief Set element to a value
		 *
		 * \param[in] iRow row number
		 * \param[in] jCol column number
		 * \param[in] input input value
		 */
		void setElem(const size_t &iRow, const size_t &jCol, const double &input);
		/** \brief Copy data from a vector to a column
		 *
		 * Copies data from a vector to a specified column. If the vector is too long, the first Nrow_ elements are used.
		 *
		 * \param[in] jCol column index (0 base)
		 * \param[in] data vector with data
		 *
		 */
		void setCol(const size_t jCol, const vector<double> data);
		/** \brief Resize matrix
		 *
		 * Resizes the matrix to the dimensions provided. All elements are set to zero.
		 *
		 * \param[in] nrow new number of rows
		 * \param[in] ncol new number of columns
		 */
		void resize(const size_t &nrow, const size_t &ncol);
		/** \brief Save matrix contents to a tab-delimited file
		 *
		 * \param[in] outFileName file name
		 */
		void save(const string &outFileName) const;

		/** \brief Vectorize the matrix
		 *
		 * Vectorize the matrix by column.
		 *
		 * \param[out] out vector of matrix elements
		 */
		void vectorize(vector<double> &out) const;
		// Linear algebra functions

		/** \brief In-place Cholesky decomposition
		 *
		 * Performs the Cholesky decomposition and stores the resulting matrix in the lower triangle of the same object.
		 */
		void chol();
		/** \brief Copy Cholesky decomposition
		 *
		 * Performs the Cholesky decomposition and stores the result in the lower triangle of the provided Matrix object. The original object is left untouched.
		 *
		 * \param[out] out object where the result is to be stored
		 */
		void chol(Matrix &out) const;
		/** \brief In-place Cholesky inverse
		 *
		 * Computes the inverse of a Cholesky decomposition and stores the resulting matrix in the same object, resulting in a symmetric matrix. The object is assumed to be a Cholesky decomposition already.
		 */
		void cholInv();
		/** \brief Copy Cholesky inverse
		 *
		 * Computes the inverse of a Cholesky decomposition and stores the result in the provided Matrix object, resulting in a symmetric matrix. The original object is left untouched. The object is assumed to be a Cholesky decomposition already.
		 *
		 * \param[out] out object where the result is to be stored
		 */
		void cholInv(Matrix &out) const;
		/** \brief Perform SVD
		 *
		 * Performs SVD and stores the \f$U\f$ vectors in a Matrix object and singular values in a C++ vector. For now, only does the _DGESVD_ from LAPACK with no \f$V^{T}\f$ matrix. The data in the object are destroyed.
		 *
		 * \param[out] U \f$U\f$ vector matrix
		 * \param[out] s singular value vector
		 */
		void svd(Matrix &U, vector<double> &s);
		/** \brief Perform "safe" SVD
		 *
		 * Performs SVD and stores the \f$U\f$ vectors in a Matrix object and singular values in a C++ vector. For now, only does the _DGESVD_ from LAPACK with no \f$V^{T}\f$ matrix. The data in the object are preserved, leading to some loss of efficiency compared to svd().
		 *
		 * \param[out] U \f$U\f$ vector matrix
		 * \param[out] s singular value vector
		 */
		void svdSafe(Matrix &U, vector<double> &s) const;
		/** \brief All eigenvalues and vectors of a symmetric matrix
		 *
		 * Interface to the _DSYEVR_ LAPACK routine. This routine is recommended as the fastest (especially for smaller matrices) in LAPACK benchmarks. It is assumed that the current object is symmetric. It is only checked for being square.
		 * The data in the relevant triangle are destroyed. If the dimensions of the output matrix and vector are smaller than necessary, they are resized. If they are larger than necessary, only the first \f$N^2\f$ and \f$N\f$ elements are used, respectively.
		 * For the matrix this means that the first \f$N\f$ columns are used if the number of rows is the same as that in the current object. Otherwise, the columns are wrapped around.
		 *
		 * \param[in] tri triangle ID ('u' for upper or 'l' for lower)
		 * \param[out] U matrix of eigenvectors
		 * \param[out] lam vector of eigenvalues in ascending order
		 *
		 */
		void eigen(const char &tri, Matrix &U, vector<double> &lam);
		/** \brief Some eigenvalues and vectors of a symmetric matrix
		 *
		 * Computes top _n_ eigenvalues and vectors of a symmetric matrix. Interface to the _DSYEVR_ LAPACK routine. This routine is recommended as the fastest (especially for smaller matrices) in LAPACK benchmarks. It is assumed that the current object is symmetric. It is only checked for being square.
		 * The data in the relevant triangle are destroyed. If the dimensions of the output matrix and vector are smaller than necessary, they are resized. If they are larger than necessary, only the first \f$N^2\f$ and \f$N\f$ elements are used, respectively.
		 * For the matrix this means that the first \f$N\f$ columns are used if the number of rows is the same as that in the current object. Otherwise, the columns are wrapped around.
		 *
		 * \param[in] tri triangle ID ('u' for upper or 'l' for lower)
		 * \param[in] n number of largest eigenvalues to compute
		 * \param[out] U matrix of eigenvectors
		 * \param[out] lam vector of eigenvalues in ascending order
		 *
		 */
		void eigen(const char &tri, const size_t &n, Matrix &U, vector<double> &lam);
		/** \brief All eigenvalues and vectors of a symmetric matrix ("safe")
		 *
		 * Interface to the _DSYEVR_ LAPACK routine. This routine is recommended as the fastest (especially for smaller matrices) in LAPACK benchmarks. It is assumed that the current object is symmetric. It is only checked for being square.
		 * The data are preserved, leading to some loss of efficiency compared to eigen(). If the dimensions of the output matrix and vector are smaller than necessary, they are resized. If they are larger than necessary, only the first \f$N^2\f$ and \f$N\f$ elements are used, respectively.
		 * For the matrix this means that the first \f$N\f$ columns are used if the number of rows is the same as that in the current object. Otherwise, the columns are wrapped around.
		 *
		 * \param[in] tri triangle ID ('u' for upper or 'l' for lower)
		 * \param[out] U matrix of eigenvectors
		 * \param[out] lam vector of eigenvalues in ascending order
		 *
		 */
		void eigenSafe(const char &tri, Matrix &U, vector<double> &lam);
		/** \brief Some eigenvalues and vectors of a symmetric matrix ("safe")
		 *
		 * Computes the top _n_ eigenvectors and values of a symmetric matrix. Interface to the _DSYEVR_ LAPACK routine. This routine is recommended as the fastest (especially for smaller matrices) in LAPACK benchmarks. It is assumed that the current object is symmetric. It is only checked for being square.
		 * The data are preserved, leading to some loss of efficiency compared to eigen(). If the dimensions of the output matrix and vector are smaller than necessary, they are resized. If they are larger than necessary, only the first \f$N^2\f$ and \f$N\f$ elements are used, respectively.
		 * For the matrix this means that the first \f$N\f$ columns are used if the number of rows is the same as that in the current object. Otherwise, the columns are wrapped around.
		 *
		 * \param[in] tri triangle ID ('u' for upper or 'l' for lower)
		 * \param[in] n number of largest eigenvalues to compute
		 * \param[out] U matrix of eigenvectors
		 * \param[out] lam vector of eigenvalues in ascending order
		 *
		 */
		void eigenSafe(const char &tri, const size_t &n, Matrix &U, vector<double> &lam);

		/** \brief In-place multiply by a design matrix from the left
		 *
		 * Performs multiplication \f$ZM\f$ where \f$Z\f$ is a design matrix (\f$Z\f$ is \f$n \times m\f$ and \f$M\f$ is \f$m \times p\f$; in practice \f$n \ge m\f$) that has one or more elements set to 1, relating the columns to rows.
		 * The function avoids actual multiplication and simply expands the \f$M\f$ matrix (in-place) to the same number of rows as \f$Z\f$. Column order in the resulting matrix is the same as in the original object.
		 *
		 * \param[in] Z \f$Z\f$ matrix
		 */
		void premultZ(const Matrix &Z);
		/** \brief Multiply by a design matrix from the left and copy result
		 *
		 * Performs multiplication \f$ZM\f$ where \f$Z\f$ is a design matrix (\f$Z\f$ is \f$n \times m\f$ and \f$M\f$ is \f$m \times p\f$; in practice \f$n \ge m\f$) that has one or more elements set to 1, relating the columns to rows.
		 * The function avoids actual multiplication and simply expands the \f$M\f$ matrix (and copies to the provided Matrix object) to the same number of rows as \f$Z\f$. Column order in the resulting matrix is the same as in the original object.
		 * If the output _Matrix_ object does not have the correct dimensions, it is resized.
		 *
		 * \param[in] Z \f$Z\f$ matrix
		 * \param[out] out output matrix
		 */
		void premultZ(const Matrix &Z, Matrix &out) const;
		/** \brief In-place multiply by the transpose of a design matrix from the left
		 *
		 * Performs multiplication \f$Z^{T}M\f$ where \f$Z\f$ is a design matrix (\f$Z\f$ is \f$m \times n\f$ and \f$M\f$ is \f$m \times p\f$; in practice \f$m \ge n\f$) that has one or more elements set to 1, relating the columns to rows. The effect is to sum the rows of \f$Y\f$ within categories represented by columns of \f$Z\f$.
		 * The function avoids actual multiplication and simply shrinks the \f$M\f$ matrix (in-place) to the same number of rows as there are columns in \f$Z\f$. Column order in the resulting matrix is the same as in the original object.
		 *
		 * \param[in] Z \f$Z\f$ matrix (not transposed)
		 */
		void premultZt(const Matrix &Z);
		/** \brief Multiply by the transpose of a design matrix from the left and copy result
		 *
		 * Performs multiplication \f$Z^{T}M\f$ where \f$Z\f$ is a design matrix (\f$Z\f$ is \f$m \times n\f$ and \f$M\f$ is \f$m \times p\f$; in practice \f$m \ge n\f$) that has one or more elements set to 1, relating the columns to rows. The effect is to sum the rows of \f$Y\f$ within categories represented by columns of \f$Z\f$.
		 * The function avoids actual multiplication and simply shrinks the \f$M\f$ matrix to the same number of rows as there are columns in \f$Z\f$. Column order in the resulting matrix is the same as in the original object.
		 * If the output _Matrix_ object does not have the correct dimensions, it is resized.
		 *
		 * \param[in] Z \f$Z\f$ matrix (not transposed)
		 * \param[out] out output matrix
		 */
		void premultZt(const Matrix &Z, Matrix &out) const;
		/** \brief In-place multiply by a design matrix from the right
		 *
		 * Performs multiplication \f$MZ\f$ where \f$Z\f$ is a design matrix (\f$Z\f$ is \f$m \times n\f$ and \f$M\f$ is \f$m \times p\f$; in practice \f$m \ge n\f$) that has one or more elements set to 1, relating the columns to rows.
		 * The function avoids actual multiplication and simply shrinks the \f$M\f$ matrix (in-place) to the number of columns the same as the number of columns in \f$Z\f$. Row order in the resulting matrix is the same as in the original object.
		 *
		 * \param[in] Z \f$Z\f$ matrix (not transposed)
		 */
		void postmultZ(const Matrix &Z);
		/** \brief Multiply by a design matrix from the right and copy result
		 *
		 * Performs multiplication \f$MZ\f$ where \f$Z\f$ is a design matrix (\f$Z\f$ is \f$m \times n\f$ and \f$M\f$ is \f$m \times p\f$; in practice \f$m \ge n\f$) that has one or more elements set to 1, relating the columns to rows.
		 * The function avoids actual multiplication and simply shrinks the \f$M\f$ matrix (in-place) to the number of columns the same as the number of columns in \f$Z\f$. Row order in the resulting matrix is the same as in the original object.
		 * If the output _Matrix_ object does not have the correct dimensions, it is resized.
		 *
		 * \param[in] Z \f$Z\f$ matrix (not transposed)
		 * \param[out] out output matrix
		 */
		void postmultZ(const Matrix &Z, Matrix &out) const;
		/** \brief In-place multiply by the transpose of a design matrix from the right
		 *
		 * Performs multiplication \f$MZ^{T}\f$ where \f$Z\f$ is a design matrix (\f$Z\f$ is \f$n \times m\f$ and \f$M\f$ is \f$m \times p\f$; in practice \f$n \ge m\f$) that has one or more elements set to 1, relating the columns to rows.
		 * The function avoids actual multiplication and simply expands the \f$M\f$ matrix (in-place) to the number of columns the same as the number of rows in \f$Z\f$. Row order in the resulting matrix is the same as in the original object.
		 *
		 * \param[in] Z \f$Z\f$ matrix (not transposed)
		 */
		void postmultZt(const Matrix &Z);
		/** \brief Multiply by the transpose of a design matrix from the right and copy result
		 *
		 * Performs multiplication \f$MZ^{T}\f$ where \f$Z\f$ is a design matrix (\f$Z\f$ is \f$n \times m\f$ and \f$M\f$ is \f$m \times p\f$; in practice \f$n \ge m\f$) that has one or more elements set to 1, relating the columns to rows.
		 * The function avoids actual multiplication and simply expands the \f$M\f$ matrix (and copies to the provided Matrix object) to the number of columns the same as the number of rows in \f$Z\f$. Row order in the resulting matrix is the same as in the original object.
		 * If the output _Matrix_ object does not have the correct dimensions, it is resized.
		 *
		 * \param[in] Z \f$Z\f$ matrix (not transposed)
		 * \param[out] out output matrix
		 */
		void postmultZt(const Matrix &Z, Matrix &out) const;

		// BLAS interface
		/** \brief Inner self crossproduct
		 *
		 * Interface for the BLAS _DSYRK_ routine. This function updates the given symmetric matrix \f$C\f$ with the operation
		 *
		 * \f$C \leftarrow \alpha A^{T}A + \beta C \f$
		 *
		 * The _char_ parameter governs which triangle of \f$C\f$ is used to store the result ('u' is upper and 'l' is lower).  If _C_ does not have the right dimensions, it is re-sized and all elements are set to zero before the operation. Otherwize, only the specified triangle is changed.
		 *
		 * \param[in] tri \f$C\f$ triangle ID
		 * \param[in] alpha the \f$\alpha\f$ parameter
		 * \param[in] beta the \f$\beta\f$ parameter
		 * \param[in,out] C the result \f$C\f$ matrix
		 */
		void syrk(const char &tri, const double &alpha, const double &beta, Matrix &C) const;
		/** \brief Outer self crossproduct
		 *
		 * Interface for the BLAS _DSYRK_ routine. This function updates the given symmetric matrix \f$C\f$ with the operation
		 *
		 * \f$C \leftarrow \alpha AA^{T} + \beta C \f$
		 *
		 * The _char_ parameter governs which triangle of \f$C\f$ is used to store the result ('u' is upper and 'l' is lower). If _C_ does not have the right dimensions, it is re-sized and all elements are set to zero before the operation. Otherwize, only the specified triangle is changed.
		 *
		 * \param[in] tri \f$C\f$ triangle ID
		 * \param[in] alpha the \f$\alpha\f$ parameter
		 * \param[in] beta the \f$\beta\f$ parameter
		 * \param[in,out] C the result \f$C\f$ matrix
		 */
		void tsyrk(const char &tri, const double &alpha, const double &beta, Matrix &C) const;
		/** \brief Multiply by symmetric matrix
		 *
		 * Multiply the _Matrix_ object by a symmetric matrix. The interface for the BLAS _DSYMM_ routine. Updates the input/output matrix \f$C\f$
		 *
		 * \f$C \leftarrow \alpha AB + \beta C \f$
		 *
		 * if _side_ is 'l' (left) and
		 *
		 * \f$C \leftarrow \alpha BA + \beta C \f$
		 *
		 * if _side_ is 'r' (right). The symmetric \f$A\f$ matrix is provided as input, the object from which the method is called is the \f$B\f$ matrix.
		 * If _C_ does not have the right dimensions, it is re-sized and all elements are set to zero before the operation. Otherwize, only the specified triangle is changed.
		 *
		 * \param[in] tri \f$A\f$ triangle ID ('u' for upper or 'l' for lower)
		 * \param[in] side multiplication side
		 * \param[in] alpha the \f$\alpha\f$ constant
		 * \param[in] symA symmetric matrix \f$A\f$
		 * \param[in] beta the \f$\beta\f$ constant
		 * \param[in,out] C the result \f$C\f$ matrix
		 *
		 */
		void symm(const char &tri, const char &side, const double &alpha, const Matrix &symA, const double &beta, Matrix &C) const;
		/** Multiply symmetric matrix by a column of another matrix
		 *
		 * Multiply the _Matrix_ object, which is symmetric, by a specified column of another matrix. An interface for the BLAS _DSYMV_ routine. Updates the input vector \f$y\f$
		 *
		 * \f$y \leftarrow \alpha AX_{\cdot j} + \beta y  \f$
		 *
		 * If the output vector is too short it is re-sized, adding zero elements as needed. If it is too long, only the first Nrow(A) elements are modified.
		 *
		 * \param[in] tri \f$A\f$ (focal object) triangle ID ('u' for upper or 'l' for lower)
		 * \param[in] alpha the \f$\alpha\f$ constant
		 * \param[in] X matrix \f$X\f$ whose column will be used
		 * \param[in] xCol column of \f$X\f$ to be used (0 base)
		 * \param[in] beta the \f$\beta\f$ constant
		 * \param[in,out] y result vector
		 *
		 */
		void symc(const char &tri, const double &alpha, const Matrix &X, const size_t &xCol, const double &beta, vector<double> &y) const;
		/** \brief General matrix multiplication
		 *
		 * Interface for the BLAS _DGEMM_ routine. Updates the input/output matrix \f$C\f$
		 *
		 * \f$ C \leftarrow \alpha op(A)op(B) + \beta C \f$
		 *
		 * where \f$op(A)\f$ is \f$A^T\f$ or \f$A\f$ if _transA_ is true or false, respectively, and similarly for \f$op(B)\f$. The object from which the method is called is \f$B\f$.
		 * If _C_ does not have the right dimensions, it is re-sized and all elements are set to zero before the operation.
		 *
		 * \param[in] transA whether \f$A\f$ should be transposed
		 * \param[in] alpha the \f$\alpha\f$ constant
		 * \param[in] A matrix \f$A\f$
		 * \param[in] transB whether \f$B\f$ should be transposed
		 * \param[in] beta the \f$\beta\f$ constant
		 * \param[in,out] C the result \f$C\f$ matrix
		 *
		 */
		void gemm(const bool &transA, const double &alpha, const Matrix &A, const bool &transB, const double &beta, Matrix &C) const;
		/** \brief Multiply a general matrix by a column of another matrix
		 *
		 * Multiply the _Matrix_ object by a specified column of another matrix. An interface for the BLAS _DGEMV_ routine. Updates the input vector \f$y\f$
		 *
		 * \f$y \leftarrow \alpha AX_{\cdot j} + \beta y  \f$
		 *
		 * or
		 *
		 * \f$y \leftarrow \alpha A^{T}X_{\cdot j} + \beta y  \f$
		 *
		 * If the output vector is too short it is re-sized, adding zero elements as needed. If it is too long, only the first Nrow(A) elements are modified.
		 *
		 * \param[in] trans whether \f$A\f$ (focal object) should be transposed
		 * \param[in] alpha the \f$\alpha\f$ constant
		 * \param[in] X matrix \f$X\f$ whose column will be used
		 * \param[in] xCol column of \f$X\f$ to be used (0 base)
		 * \param[in] beta the \f$\beta\f$ constant
		 * \param[in,out] y result vector
		 *
		 */
		void gemc(const bool &trans, const double &alpha, const Matrix &X, const size_t &xCol, const double &beta, vector<double> &y) const;

		// Sampling functions
		/** \brief Shuffle columns
		 *
		 * Suffle the columns of the current object and return a matrix with of the same size but column order randomly permuted.
		 *
		 * \return permuted `Matrix` object
		 */
		Matrix colShuffle() const;
		/** \brief Shuffle rows
		 *
		 * Suffle the rows of the current object and return a matrix with of the same size but row order randomly permuted.
		 *
		 * \return permuted `Matrix` object
		 */
		Matrix rowShuffle() const;

		// Overloaded operators
		/** \brief Hadamard matrix product
		 *
		 * \param[in] m right matrix
		 * \return Matrix result
		 *
		 */
		Matrix operator*(const Matrix &m) const;
		/** \brief Matrix-scalar product
		 *
		 * \param[in] scal scalar
		 * \return Matrix result
		 *
		 */
		Matrix operator*(const double &scal) const;
		/** \brief Entrywise matrix division
		 *
		 * \param[in] m right matrix
		 * \return Matrix result
		 *
		 */
		Matrix operator/(const Matrix &m) const;
		/** \brief Matrix-scalar division
		 *
		 * \param[in] scal scalar
		 * \return Matrix result
		 *
		 */
		Matrix operator/(const double &scal) const;
		/** \brief Entrywise matrix addition
		 *
		 * \param[in] m right matrix
		 * \return Matrix result
		 *
		 */
		Matrix operator+(const Matrix &m) const;
		/** \brief Matrix-scalar addition
		 *
		 * \param[in] scal scalar
		 * \return Matrix result
		 *
		 */
		Matrix operator+(const double &scal) const;
		/** \brief Entrywise matrix subtraction
		 *
		 * \param[in] m right matrix
		 * \return Matrix result
		 *
		 */
		Matrix operator-(const Matrix &m) const;
		/** \brief Matrix-scalar subtraction
		 *
		 * \param[in] scal scalar
		 * \return Matrix result
		 *
		 */
		Matrix operator-(const double &scal) const;

		// compound assignment operators
		/** \brief Matrix-scalar compound addition
		 *
		 * \param[in] scal scalar
		 * \return Matrix result
		 *
		 */
		Matrix& operator+=(const double &scal);
		/** \brief Matrix-scalar compound product
		 *
		 * \param[in] scal scalar
		 * \return Matrix result
		 *
		 */
		Matrix& operator*=(const double &scal);
		/** \brief Matrix-scalar compound subtraction
		 *
		 * \param[in] scal scalar
		 * \return Matrix result
		 *
		 */
		Matrix& operator-=(const double &scal);
		/** \brief Matrix-scalar compound division
		 *
		 * \param[in] scal scalar
		 * \return Matrix result
		 *
		 */
		Matrix& operator/=(const double &scal);

		// column- and row-wise operations
		/** \brief Row means
		 *
		 * Calculates row means and stores them in the provided vector. If vector length is smaller than necessary, the vector is expanded. Otherwise, the first \f$N_{row}\f$ elements are used.
		 *
		 * \param[out] means vector of means
		 *
		 */
		void rowMeans(vector<double> &means) const;
		/** \brief Column means
		 *
		 * Calculates column means and stores them in the provided vector. If vector length is smaller than necessary, the vector is expanded. Otherwise, the first \f$N_{col}\f$ elements are used.
		 *
		 * \param[out] means vector of means
		 *
		 */
		void colMeans(vector<double> &means) const;
		/** \brief Row sums
		 *
		 * Calculates row sums and stores them in the provided vector. If vector length is smaller than necessary, the vector is expanded. Otherwise, the first \f$N_{row}\f$ elements are used.
		 *
		 * \param[out] sums vector of sums
		 *
		 */
		void rowSums(vector<double> &sums) const;
		/** \brief Column sums
		 *
		 * Calculates column sums and stores them in the provided vector. If vector length is smaller than necessary, the vector is expanded. Otherwise, the first \f$N_{col}\f$ elements are used.
		 *
		 * \param[out] sums vector of sums
		 *
		 */
		void colSums(vector<double> &sums) const;
		/** \brief Multiply rows by a vector
		 *
		 * Entry-wise multiplication of each row by the provided vector. The current object is modified.
		 *
		 * \param[in] scalars vector of scalars to use for multiplication
		 *
		 */
		void rowMultiply(const vector<double> &scalars);
		/** \brief Multiply a row by a scalar
		 *
		 * Entry-wise multiplication of a given row by the provided scalar. The current object is modified.
		 *
		 * \param[in] scalar scalar to use for multiplication
		 * \param[in] iRow row index
		 *
		 */
		void rowMultiply(const double &scalar, const size_t &iRow);
		/** \brief Multiply columns by a vector
		 *
		 * Entry-wise multiplication of each column by the provided vector. The current object is modified.
		 *
		 * \param[in] scalars vector of scalars to use for multiplication
		 *
		 */
		void colMultiply(const vector<double> &scalars);
		/** \brief Multiply a column by a scalar
		 *
		 * Entry-wise multiplication of a given column by the provided scalar. The current object is modified.
		 *
		 * \param[in] scalar scalar to use for multiplication
		 * \param[in] jCol column index
		 *
		 */
		void colMultiply(const double &scalar, const size_t &jCol);
		/** \brief Divide rows by a vector
		 *
		 * Entry-wise division of each row by the provided vector. The current object is modified.
		 *
		 * \param[in] scalars vector of scalars to use for division
		 *
		 */
		void rowDivide(const vector<double> &scalars);
		/** \brief Divide a row by a scalar
		 *
		 * Entry-wise division of a given row by the provided scalar. The current object is modified.
		 *
		 * \param[in] scalar scalar to use for division
		 * \param[in] iRow row index
		 *
		 */
		void rowDivide(const double &scalar, const size_t &iRow);
		/** \brief Divide columns by a vector
		 *
		 * Entry-wise division of each column by the provided vector. The current object is modified.
		 *
		 * \param[in] scalars vector of scalars to use for division
		 *
		 */
		void colDivide(const vector<double> &scalars);
		/** \brief Divide a column by a scalar
		 *
		 * Entry-wise division of a given column by the provided scalar. The current object is modified.
		 *
		 * \param[in] scalar scalar to use for division
		 * \param[in] jCol column index
		 *
		 */
		void colDivide(const double &scalar, const size_t &jCol);
		/** \brief Add a vector to rows
		 *
		 * Entry-wise addition of a vector to each row. The current object is modified.
		 *
		 * \param[in] scalars vector of scalars to use for addition
		 *
		 */
		void rowAdd(const vector<double> &scalars);
		/** \brief Add a scalar to a row
		 *
		 * Entry-wise addition of a scalar to the given row. The current object is modified.
		 *
		 * \param[in] scalar scalar to use for addition
		 * \param[in] iRow row index
		 *
		 */
		void rowAdd(const double &scalar, const size_t &iRow);
		/** \brief Add a vector to columns
		 *
		 * Entry-wise addition of a vector to each column. The current object is modified.
		 *
		 * \param[in] scalars vector of scalars to use for addition
		 *
		 */
		void colAdd(const vector<double> &scalars);
		/** \brief Add a scalar to a column
		 *
		 * Entry-wise addition of a scalar to the given column. The current object is modified.
		 *
		 * \param[in] scalar scalar to use for addition
		 * \param[in] jCol column index
		 *
		 */
		void colAdd(const double &scalar, const size_t &jCol);
		/** \brief Subtract a vector from rows
		 *
		 * Entry-wise subtraction of a vector from each row. The current object is modified.
		 *
		 * \param[in] scalars vector of scalars to use for subtraction
		 *
		 */
		void rowSub(const vector<double> &scalars);
		/** \brief Subtract a scalar from a row
		 *
		 * Entry-wise subtraction of a scalar from the given row. The current object is modified.
		 *
		 * \param[in] scalar scalar to use for subtraction
		 * \param[in] iRow row index
		 *
		 */
		void rowSub(const double &scalar, const size_t &iRow);
		/** \brief Subtract a vector from columns
		 *
		 * Entry-wise subtraction of a vector from each column. The current object is modified.
		 *
		 * \param[in] scalars vector of scalars to use for subtraction
		 *
		 */
		void colSub(const vector<double> &scalars);
		/** \brief Subtract a scalar from a column
		 *
		 * Entry-wise subtraction of a scalar from the given column. The current object is modified.
		 *
		 * \param[in] scalar scalar to use for subtraction
		 * \param[in] jCol column index
		 *
		 */
		void colSub(const double &scalar, const size_t &jCol);
		/** \brief Append columns of a matrix
		 *
		 * Columns of a provided matrix are appended after the last column of the current object. The object is expanded accordingly. Number of rows does not change.
		 *
		 * \param[in] cols a _Matrix_ object with columns to append
		 *
		 */
		void appendCol(const Matrix &cols);
		/** \brief Append rows of a matrix
		 *
		 * Rows of a provided matrix are appended after the last row of the current object. The object is expanded accordingly. Number of columns does not change.
		 *
		 * \param[in] rows a _Matrix_ object with rows to append
		 *
		 */
		void appendRow(const Matrix &rows);
		/** \brief Drop left columns
		 *
		 * Delete columns to the left of the one marked by the provided index.
		 *
		 * \param[in] newFirst index of the new first column
		 *
		 */
		void dropLeftCols(const size_t &newFirst);
		/** \brief Drop right columns
		 *
		 * Delete columns to the right of the one marked by the provided index.
		 *
		 * \param[in] newLast index of the new last column
		 *
		 */
		void dropRightCols(const size_t &newLast);
		/** \brief Drop top rows
		 *
		 * Delete rows higher than the one marked by the provided index.
		 *
		 * \param[in] newTop index of the new top row
		 *
		 */
		void dropTopRows(const size_t &newTop);
		/** \brief Drop bottom rows
		 *
		 * Delete rows lower than the one marked by the provided index.
		 *
		 * \param[in] newBottom index of the new bottom row
		 *
		 */
		void dropBottomRows(const size_t &newBottom);
	private:
		double *data_;
		size_t Nrow_;
		size_t Ncol_;
	};

}


#endif /* locMatrix_hpp */
