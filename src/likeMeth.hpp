//
//  likeMeth.hpp
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
 * This is the project header file containing class definitions and interface documentation.
 *
 */

#ifndef likeMeth_hpp
#define likeMeth_hpp

#include <vector>

#include "locMatrix.hpp"

using std::vector;
using locMatrix::Matrix;

class EmmREML;
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
private:
	/// \f$\eta^2\f$ from Kang _et al_ equation (7). Points to a Matrix object with each trait as a column.
	const Matrix *_etaSq;
	/// \f$\lambda\f$ from Kang _et al_ equation (7). Points to a vector of eigenvalues.
	const vector<double> *_lambda;
	/// column ID of the etaSq matrix
	size_t _jCol;

public:

	/// Default constructor
	EmmREML() : _etaSq(nullptr), _lambda(nullptr), _jCol(0) {};
	/** \brief Constructor
	 *
	 * Sets up the object.
	 *
	 * \param[in] etaSq address of a Matrix object with \f$\eta^2\f$
	 * \param[in] lambda address of a Matrix object with \f$\lambda\f$
	 * \param[in] jCol phenotype column index
	 *
	 */
	EmmREML(const Matrix &etaSq, const vector<double> &lambda, const size_t &jCol) : _etaSq(&etaSq), _lambda(&lambda), _jCol(jCol) {};

	/// Destructor
	~EmmREML(){_etaSq = nullptr; _lambda = nullptr; };

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
	void setColID(const size_t &j) {_jCol = j; };

};

/** \brief A SNP block functor class
 *
 * A class to facilitate GWA multithreading. Points to a block of SNPs and runs _snpReg()_ on each locus.
 *
 */
class SNPblock {
private:
	/// Pointer to the response matrix
	const Matrix *_rsp;
	/// Pointer to the SNP array
	const int *_snp;
	/// Pointer to the \f$-\log_{10}p\f$ array (the results)
	double *_lPval;
	/// Number of SNPs in a block
	size_t _blockSize;
	/// Start position of the block in the SNP array
	size_t _blockStart;

public:
	/// Default constructor
	SNPblock() : _rsp(nullptr), _snp(nullptr), _lPval(nullptr), _blockSize(0) {};
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
	SNPblock(const Matrix &rsp, const int *snpArr, const size_t &bStart, const size_t &bSize, double *lpvArr): _rsp(&rsp), _snp(snpArr), _blockSize(bSize), _blockStart(bStart), _lPval(lpvArr) {};

	/// Destructor
	~SNPblock(){ _rsp = nullptr; _snp = nullptr; _lPval = nullptr; };

	/// Copy constructor (not implemented)
	SNPblock(const SNPblock &) = delete;
	/// Copy assignment operator (not implemented)
	SNPblock & operator=(const SNPblock &) = delete;
	/// Move constructor
	SNPblock(SNPblock &&in) : _rsp(in._rsp), _snp(in._snp), _lPval(in._lPval), _blockSize(in._blockSize), _blockStart(in._blockStart) { in._rsp = nullptr; in._snp = nullptr; in._lPval = nullptr; };

	/** \brief Function operator
	 *
	 * Performs GWA on each SNP in the block.
	 *
	 */
	void operator()();
};


#endif /* likeMeth_hpp */
