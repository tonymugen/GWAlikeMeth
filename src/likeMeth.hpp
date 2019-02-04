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
 * This is the project header file containing class definitions and interface documentation.
 *
 */

#ifndef likeMeth_hpp
#define likeMeth_hpp

#include <vector>

#include "locMatrix.hpp"

using std::vector;
using std::move;
using locMatrix::Matrix;

class EmmREML;
class MixedModel;
class SNPblock;

/** \brief Solve a mixed model with no fixed effect covariates
 *
 * Solves a mixed model
 *
 * \f$ Y = \mu + Zu + e \f$
 *
 * \f$ u \sim \textrm{N}\left(0, \sigma^2_gK \right)\f$
 *
 * for each trait in the trait matrix \f$Y\f$ separately, but calculating common parameters only once.
 *
 * \param[in] Y phenotype data
 * \param[in] K relationship matrix
 * \param[in] Z design matrix
 * \param[out] u  genome-estimated breeding values; input dimensions can be arbitrary: the matrix will be re-sized
 * \param[out] mu vector of intercepts
 *
 */
void solveMM(const Matrix &Y, const Matrix &K, const Matrix &Z, Matrix &u, vector<double> &mu);
/** \brief Solve a mixed model with fixed effect covariates
 *
 * Solves a mixed model
 *
 * \f$ Y = X\beta + Zu + e \f$
 *
 * \f$ u \sim \textrm{N}\left(0, \sigma^2_gK \right)\f$
 *
 * for each trait in the trait matrix \f$Y\f$ separately, but calculating common parameters only once. The matrix \f$X\f$ contains the fixed-effect covariates and the intercept. The intercept column is added to the provided matrix with the function, so it should not be included in the parameter.
 * \f$X\f$ can have the number of rows equal to either the number of lines (rows in \f$K\f$ or columns in \f$Z\f$) or the number of data points (rows of \f$Y\f$ or \f$Z\f$). Note also that the \f$X^{T}X\f$ matrix has to be non-singular. Simple dimension checks are made to eleminate obvious singularity, but the user has to make sure \f$X\f$ is full-rank.
 *
 * \param[in] Y phenotype data
 * \param[in] K relationship matrix
 * \param[in] Z design matrix
 * \param[in] X fixed effect matrix
 * \param[out] u  genome-estimated breeding values; input dimensions can be arbitrary: the matrix will be re-sized
 * \param[out] beta matrix of fixed effects
 *
 */
void solveMM(const Matrix &Y, const Matrix &K, const Matrix &Z, const Matrix &X, Matrix &u, Matrix &beta);


/** \brief SNP regression
 *
 * Estimates regression \f$ -\log_{10}p \f$ for a SNP with missing genotype data and multiple traits in a table. Genotypes should be coded as (0,1,2) for homozygous, hetereozygous and homozygous alternative. It will run faster of the major allele homozygotes are coded as 0.
 * Each trait is treated separately but it helps to include multiple traits because some caclulations are common and can performed only once.
 *
 * \param[in] Y response matrix with traits as columns and lines (represented once each) as rows
 * \param[in] snp array of genotypes
 * \param[in] misTok missing value identifier; should be distinct from (0,1,2)
 * \param[out] lPval \f$-\log_{10}p\f$
 *
 */
void snpReg(const Matrix &Y, const int *snp, const int &misTok, double *lPval);

/** \brief EMMA REML functor class
 *
 * Class that implements the REML function from Kang _et al_ (2008). Used in the likelihood maximization step.
 *
 */
class EmmREML {
public:

	/// Default constructor
	EmmREML() : etaSq_(nullptr), lambda_(nullptr), jCol_(0) {};
	/** \brief Constructor
	 *
	 * Sets up the object.
	 *
	 * \param[in] etaSq address of a Matrix object with \f$\eta^2\f$
	 * \param[in] lambda address of a Matrix object with \f$\lambda\f$
	 * \param[in] jCol phenotype column index
	 *
	 */
	EmmREML(const Matrix &etaSq, const vector<double> &lambda, const size_t &jCol) : etaSq_(&etaSq), lambda_(&lambda), jCol_(jCol) {};

	/// Destructor
	~EmmREML(){etaSq_ = nullptr; lambda_ = nullptr; };

	/// Copy constructor (not implemented)
	EmmREML(const EmmREML &) = delete;

	/// Assignement operator (not implemented)
	EmmREML & operator=(const EmmREML &) = delete;

	/** \brief Function operator
	 *
	 * Does the restricted likelihood calculation for each value of \f$\delta = \frac{\sigma^2_e}{\sigma^2_g}\f$ and the given columns of the \f$\eta^2\f$ and \f$\lambda\f$ matrices.
	 *
	 * \param[in] delta \f$\delta\f$ value
	 *
	 */
	double operator()(const double &delta);

	/** \brief Set column index
	 *
	 *\param[in] j column index
	 */
	void setColID(const size_t &j) {jCol_ = j; };
private:
	/// \f$\eta^2\f$ from Kang _et al_ equation (7). Points to a Matrix object with each trait as a column.
	const Matrix *etaSq_;
	/// \f$\lambda\f$ from Kang _et al_ equation (7). Points to a vector of eigenvalues.
	const vector<double> *lambda_;
	/// column ID of the etaSq matrix
	size_t jCol_;


};

/** \brief Mixed model
 *
 * Constructors solve a mixed model given inputs. Public functions do GWA.
 *
 */
class MixedModel {
	public:
		/** \brief Default constructor  */
		MixedModel() : Y_{Matrix()}, K_{Matrix()}, Z_{Matrix()}, X_{Matrix()}, u_{Matrix()}, beta_{Matrix()}, delta_{0.0} {};
		/** \brief Constructor with no fixed effect
		 *
		 * \param[in] yvec vectorized (by column) response matrix
		 * \param[in] kvec vectorized relationship matrix
		 * \param[in] repFac factor relating genotypes to replicates
		 * \param[in] d number of traits
		 * \param[in] Ngen number of genotypes
		 * \param[in] N number of data points
		 *
		 */
		MixedModel(const vector<double> &yvec, const vector<double> &kvec, const vector<double> &repFac, const size_t &d, const size_t &Ngen, const size_t &N);
		/** \brief Constructor including a fixed effect
		 *
		 * \param[in] yvec vectorized (by column) response matrix
		 * \param[in] kvec vectorized relationship matrix
		 * \param[in] repFac factor relating genotypes to replicates
		 * \param[in] xvec vectorized (by column) matrix of fixed effects
		 * \param[in] d number of traits
		 * \param[in] Ngen number of genotypes
		 * \param[in] N number of data points
		 *
		 */
		MixedModel(const vector<double> &yvec, const vector<double> &kvec, const vector<double> &repFac, const vector<double> &xvec, const size_t &d, const size_t &Ngen, const size_t &N);

		/// Destructor
		~MixedModel(){};

		/// Copy constructor (not implemented)
		MixedModel(const MixedModel &in) = delete;
		/// Copy assignement (not implemented)
		MixedModel & operator=(const MixedModel &in) = delete;
		/** \brief Move constructor
		 *
		 * \param[in] in object to be moved
		 */
		MixedModel(MixedModel &&in) : Y_{move(in.Y_)}, K_{move(in.K_)}, Z_{move(in.Z_)}, X_{move(in.X_)}, u_{move(in.u_)}, delta_{in.delta_} {};
		/** \brief Move assignment operator (not impplemented) */
		MixedModel & operator=(MixedModel &&in) = delete;

		/** \brief Access the random effects
		 *
		 * \param[out] u random effect matrix
		 */
		void ranef(Matrix &u) const;
		/** \brief Access the fixed effects
		 *
		 * Returns empty matrix if there are no fixed effects.
		 *
		 * \param[out] beta fixed effect matrix
		 */
		void fixef(Matrix &beta) const;
		/** \brief Marker heritability
		 *
		 * \return marker heritability \f$ h^2_\textsc{m} \f$
		 */
		double hSq() const { return 1.0/(1.0 + delta_); };
	private:
		/** \brief Responce matrix */
		const Matrix Y_;
		/** \brief Relationship matrix */
		const Matrix K_;
		/** \brief Design matrix */
		const Matrix Z_;
		/** brief Fixed effects matrix */
		const Matrix X_;
		/** \brief BLUP matrix */
		Matrix u_;
		/** \brief Fixed effect matrix */
		Matrix beta_;
		/** \brief Delta value
		 *
		 * \f$ \delta = \dfrac{\sigma^2_e}{\sigma^2_g}\f$
		 *
		 */
		double delta_;

};

/** \brief A SNP block functor class
 *
 * A class to facilitate GWA multithreading. Points to a block of SNPs and runs snp_Reg()_ on each locus.
 *
 */
class SNPblock {
public:
	/// Default constructor
	SNPblock() : rsp_{nullptr}, snp_{nullptr}, lPval_{nullptr}, blockSize_{0} {};
	/** \brief Constructor
	 *
	 * Sets up the pointers to the data and results objects
	 *
	 * \param[in] rsp address of the response matrix
	 * \param[in] snpArr pointer to the SNP array
	 * \param[in] bStart block start position
	 * \param[in] bSize block size (number of SNPs)
	 * \param[out] lpvArr pointer to the array where the \f$-\log_{10}p\f$ are stored
	 *
	 */
	SNPblock(const Matrix &rsp, const int *snpArr, const size_t &bStart, const size_t &bSize, double *lpvArr): rsp_(&rsp), snp_(snpArr), blockSize_(bSize), blockStart_(bStart), lPval_(lpvArr) {};

	/// Destructor
	~SNPblock(){ rsp_ = nullptr; snp_ = nullptr; lPval_ = nullptr; };

	/// Copy constructor (not implemented)
	SNPblock(const SNPblock &) = delete;
	/// Copy assignment operator (not implemented)
	SNPblock & operator=(const SNPblock &) = delete;
	/// Move constructor
	SNPblock(SNPblock &&in) : rsp_{move(in.rsp_)}, snp_{move(in.snp_)}, lPval_{move(in.lPval_)}, blockSize_{move(in.blockSize_)}, blockStart_{move(in.blockStart_)} { in.rsp_ = nullptr; in.snp_ = nullptr; in.lPval_ = nullptr; };

	/** \brief Function operator
	 *
	 * Performs GWA on each SNP in the block.
	 *
	 */
	void operator()();
private:
	/// Pointer to the response matrix
	const Matrix *rsp_;
	/// Pointer to the SNP array
	const int *snp_;
	/// Pointer to the \f$-\log_{10}p\f$ array (the results)
	double *lPval_;
	/// Number of SNPs in a block
	size_t blockSize_;
	/// Start position of the block in the SNP array
	size_t blockStart_;

};


#endif /* likeMeth_hpp */
