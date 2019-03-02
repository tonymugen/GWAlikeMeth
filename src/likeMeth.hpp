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
using BayesicSpace::Matrix;



namespace BayesicSpace {
	class EmmREML;
	class MixedModel;
	class SNPblock;

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
		friend class SNPblock;

		public:
			/** \brief Default constructor  */
			MixedModel() : misTok_{0} {};
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
			MixedModel(const vector<double> &yvec, const vector<double> &kvec, const vector<size_t> &repFac, const size_t &d, const size_t &Ngen, const size_t &N);
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
			MixedModel(const vector<double> &yvec, const vector<double> &kvec, const vector<size_t> &repFac, const vector<double> &xvec, const size_t &d, const size_t &Ngen, const size_t &N);
			/** \brief Constructor with SNPs but no fixed effect
			 *
			 * \param[in] yvec vectorized (by column) response matrix
			 * \param[in] kvec vectorized relationship matrix
			 * \param[in] repFac factor relating genotypes to replicates
			 * \param[in] d number of traits
			 * \param[in] Ngen number of genotypes
			 * \param[in] N number of data points
			 * \param[in] snps vectorized (by column) SNP matrix
			 * \param[in] misTok token representing missing genotype data
			 * \param[out] lPval address of the voectorized \f$ -\log_{10}p \f$ matrix
			 *
			 */
			MixedModel(const vector<double> &yvec, const vector<double> &kvec, const vector<size_t> &repFac, const size_t &d, const size_t &Ngen, const size_t &N, const vector<int32_t> *snps, const int32_t &misTok, vector<double> *lPval) : MixedModel(yvec, kvec, repFac, d, Ngen, N) {snps_ = snps; misTok_ = misTok; lPval_ = lPval; };
			/** \brief Constructor including SNPs and a fixed effect
			 *
			 * \param[in] yvec vectorized (by column) response matrix
			 * \param[in] kvec vectorized relationship matrix
			 * \param[in] repFac factor relating genotypes to replicates
			 * \param[in] xvec vectorized (by column) matrix of fixed effects
			 * \param[in] d number of traits
			 * \param[in] Ngen number of genotypes
			 * \param[in] N number of data points
			 * \param[in] snps vectorized (by column) SNP matrix
			 * \param[in] misTok token representing missing genotype data
			 * \param[out] lPval address of the voectorized \f$ -\log_{10}p \f$ matrix
			 *
			 */
			MixedModel(const vector<double> &yvec, const vector<double> &kvec, const vector<size_t> &repFac, const vector<double> &xvec, const size_t &d, const size_t &Ngen, const size_t &N, const vector<int32_t> *snps, const int32_t &misTok, vector<double> *lPval) : MixedModel(yvec, kvec, repFac, xvec, d, Ngen, N) {snps_ = snps; misTok_ = misTok; lPval_ = lPval; };
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
			MixedModel(MixedModel &&in) : Y_{move(in.Y_)}, K_{move(in.K_)}, Z_{move(in.Z_)}, X_{move(in.X_)}, u_{move(in.u_)}, delta_{move(in.delta_)}, snps_{move(in.snps_)}, misTok_{in.misTok_}, lPval_{move(in.lPval_)} {in.snps_ = nullptr; in.lPval_ = nullptr; };
			/** \brief Move assignment operator (not impplemented) */
			MixedModel & operator=(MixedModel &&in) = delete;

			/** \brief Access the random effects
			 *
			 * \param[out] u random effect matrix
			 */
			void ranef(Matrix &u) const {u = u_; };
			/** \brief Access the fixed effects
			 *
			 * Returns empty matrix if there are no fixed effects.
			 *
			 * \param[out] beta fixed effect matrix
			 */
			void fixef(Matrix &beta) const {beta = beta_; };
			/** \brief Marker heritability
			 *
			 * \param[out] marker heritability \f$ h^2_\textsc{m} \f$
			 */
			void hSq(vector<double> &out) const;
			/** \brief Genome-wide association
			 *
			 * Estimates regression \f$ -\log_{10}p \f$ for SNPs with missing genotype data and multiple traits in a table. Genotypes should be coded as (0,1,2) for homozygous, hetereozygous and homozygous alternative. It should run faster of the major allele homozygotes are coded as 0.
			 * Each trait is treated separately but it helps to include multiple traits because some calculations are common and can be performed only once. The \f$ -\log_{10}p \f$ matrix has each trait in a column.
			 *
			 */
			void gwa();
			/** \brief Genome-wide association with FDR
			 *
			 * The same as `gwa()`, but does permutations to calculate emprical FDR for each SNP.
			 *
			 * \param[in] nPer number of permutations
			 * \param[out] fdr vectorized (by column) matrix of FDR values
			 *
			 */
			void gwa(const uint32_t &nPer, vector<double> &fdr);
		private:
			/** \brief Responce matrix */
			Matrix Y_;
			/** \brief Relationship matrix */
			const Matrix K_;
			/** \brief Design matrix */
			Matrix Z_;
			/** brief Fixed effects matrix */
			const Matrix X_;
			/** \brief BLUP matrix */
			Matrix u_;
			/** \brief Fixed effect matrix */
			Matrix beta_;
			/** \brief Delta values
			 *
			 * \f$ \delta = \dfrac{\sigma^2_e}{\sigma^2_g}\f$, one for each trait.
			 *
			 */
			vector<double> delta_;
			/** \brief Pointer to a SNP matrix
			 *
			 * Pointer to a vectorized (by column) SNP matrix. Individuals are in rows, SNPs in columns.
			 */
			const vector<int32_t> *snps_;
			/** \brief Missing genotype token */
			int32_t misTok_;
			/** \brief Pointer to the \f$ -\log_{10}p \f$ matrix
			 *
			 * Vectorized (by column; traits are in columns) matrix.
			 */
			vector<double> *lPval_;
			/** \brief Single SNP regression
			 *
			 * \param[in] idx index of the SNP
			 * \param[in] Nsnp number of SNPs
			 *
			 */
			void oneSNP_(const size_t &idx, const size_t &Nsnp);
			/** \brief Single SNP regression for permutations
			 *
			 * The permuted \f$ -\log_{10}p \f$ vector has the values for each trait continuous in memory, each permutation after the other: trait1:per1|per2|...|perN->trait2:per1|per2|...|perN->...
			 *
			 * \param[in] snpIdx index of the current SNP
			 * \param[in] perOff permutation offset (`nPer*Nsnp`)
			 * \param[in] snpOff SNP offset for each trait (`Nsnp*perIdx`)
			 * \param[out] lPval vector of permutation \f$ -\log_{10}p \f$
			 *
			 */
			void oneSNP_(const size_t &snpIdx, const size_t &perOff, const size_t &snpOff, vector<double> &lPval);
	};

	/** \brief A SNP block functor class
	 *
	 * A class to facilitate GWA multithreading in the `MixedModel` class. Points to a block of SNPs and runs `oneSNP_()` on each locus.
	 *
	 */
	class SNPblock {
	public:
		/// Default constructor
		SNPblock() : mmObj_{nullptr}, blockSize_{0}, blockStart_{0}, plPval_{nullptr}, perOff_{0}, snpOff_{0} {};
		/** \brief Constructor
		 *
		 * Sets up the pointer to the `MixedModel` object calling the current instance of this functor.
		 *
		 * \param[in,out] parent address of the `MixedModel` object calling this functor
		 * \param[in] bStart block start position
		 * \param[in] bSize block size (number of SNPs)
		 *
		 */
		SNPblock(MixedModel &parent, const size_t &bStart, const size_t &bSize): mmObj_{&parent}, blockSize_{bSize}, blockStart_{bStart}, plPval_{nullptr}, perOff_{0}, snpOff_{0} {};
		/** \brief Constructor for permutations
		 *
		 * Sets up the pointer to the `MixedModel` object calling the current instance of this functor, adding the info to work on permuted data.
		 *
		 * \param[in,out] parent address of the `MixedModel` object calling this functor
		 * \param[in] bStart block start position
		 * \param[in] bSize block size (number of SNPs)
		 * \param[out] plPval address of the permuted \f$-\log_{10}p \f$ vector
		 * \param[in] perOff permutation offset (`Nsnp*nPer`)
		 * \param[in] snpOff SNP offset (`Nsnp*perIdx`)
		 *
		 */
		SNPblock(MixedModel &parent, const size_t &bStart, const size_t &bSize, vector<double> &plPval, const size_t &perOff, const size_t &snpOff): mmObj_{&parent}, blockSize_{bSize}, blockStart_{bStart}, plPval_{&plPval}, perOff_{perOff}, snpOff_{snpOff} {};

		/// Destructor
		~SNPblock(){ mmObj_ = nullptr; };

		/// Copy constructor (not implemented)
		SNPblock(const SNPblock &) = delete;
		/// Copy assignment operator (not implemented)
		SNPblock & operator=(const SNPblock &) = delete;
		/// Move constructor
		SNPblock(SNPblock &&in) : mmObj_{move(in.mmObj_)}, blockSize_{move(in.blockSize_)}, blockStart_{move(in.blockStart_)}, plPval_{move(in.plPval_)}, perOff_{move(in.perOff_)}, snpOff_{move(in.snpOff_)} { in.mmObj_ = nullptr; in.plPval_ = nullptr; };

		/** \brief Function operator
		 *
		 * Performs GWA on each SNP in the block.
		 *
		 */
		void operator()();
	private:
		/// Pointer to a `MixedModel` object
		MixedModel *mmObj_;
		/// Number of SNPs in a block
		size_t blockSize_;
		/// Start position of the block in the SNP array
		size_t blockStart_;
		/// Pointer to permuted \f$ -\log_{10}p \f$ vector
		vector<double> *plPval_;
		/** \brief Permutation offset
		 *
		 * The `Nsnp*nPer` value fed to `oneSNP()` from the `MixedModel` class.
		 */
		size_t perOff_;
		/** \brief SNP offset
		 *
		 * The `Nsnp*perIdx` value fed to `oneSNP()` from the `MixedModel` class.
		 */
		size_t snpOff_;
	};

}
#endif /* likeMeth_hpp */

