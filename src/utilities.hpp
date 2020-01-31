//
//  utilities.hpp
//
//
//  Created by Tony Greenberg on 12/21/16.
//  Copyright © 2016 Tony Greenberg. All rights reserved.
//

/// Miscellaneous functions and algorithms
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2016 Anthony J. Greenberg
 * \version 0.1
 *
 * This is the project header file containing function definitions and constants.
 *
 */

#ifndef utilities_hpp
#define utilities_hpp

#include <vector>
#include <string>
#include <utility>
#include <limits>
#include <cmath>
#include <cstdint> // for uint64_t
#include <stack>

using std::vector;
using std::string;
using std::swap;
using std::numeric_limits;
using std::signbit;
using std::stack;

namespace BayesicSpace {
	/// The definition of \f$\pi\f$
	const double BS_PI = 3.14159265358979323846264338328;
	/// Machine \f$\epsilon\f$
	const double BS_EPS   = numeric_limits<double>::epsilon();
	/// Tiny value to guard agains underflow
	const double BS_FPMIN = numeric_limits<double>::min()/BS_EPS;

	/** \brief Swap two `uint64_t` values
	 *
	 * Uses the three XORs trick to swap two integers. Safe if the variables happen to refer to the same address.
	 *
	 * \param[in,out] i first integer
	 * \param[in,out] j second integer
	 */
	void swapXOR(uint64_t &i, uint64_t &j){
		if (&i != &j) { // no move needed if this is actually the same variable
			i ^= j;
			j ^= i;
			i ^= j;
		}
	}
	/** \brief Swap two `in64_t` values
	 *
	 * Uses the three XORs trick to swap two integers. Safe if the variables happen to refer to the same address.
	 *
	 * \param[in,out] i first integer
	 * \param[in,out] j second integer
	 */
	void swapXOR(int64_t &i, int64_t &j){
		if (&i != &j) { // no move needed if this is actually the same variable
			i ^= j;
			j ^= i;
			i ^= j;
		}
	}
	/** \brief Swap two `uint32_t` values
	 *
	 * Uses the three XORs trick to swap two integers. Safe if the variables happen to refer to the same address.
	 *
	 * \param[in,out] i first integer
	 * \param[in,out] j second integer
	 */
	void swapXOR(uint32_t &i, uint32_t &j){
		if (&i != &j) { // no move needed if this is actually the same variable
			i ^= j;
			j ^= i;
			i ^= j;
		}
	}
	/** \brief Swap two `int32_t` values
	 *
	 * Uses the three XORs trick to swap two integers. Safe if the variables happen to refer to the same address.
	 *
	 * \param[in,out] i first integer
	 * \param[in,out] j second integer
	 */
	void swapXOR(int32_t &i, int32_t &j){
		if (&i != &j) { // no move needed if this is actually the same variable
			i ^= j;
			j ^= i;
			i ^= j;
		}
	}
	/** \brief Mean of a C array
	 *
	 * Calculates a mean of an array.  Uses the recursive algorithm for numerical stability.
	 *
	 * \param[in] arr array to average
	 * \param[in] len array length
	 *
	 */
	double mean(const double *arr, const size_t &len){
		double mean = 0.0;

		for (size_t i = 0; i < len; i++) {
			mean += (arr[i] - mean)/static_cast<double>(i + 1);
		}
		return mean;
	}

	/** \brief Mean of a C array with stride
	 *
	 * Calculates a mean of an array with stride (i.e., using every _stride_ element).  Uses the recursive algorithm for numerical stability.
	 *
	 * \param[in] arr array to average
	 * \param[in] len array length
	 * \param[in] stride stride
	 *
	 */
	double mean(const double *arr, const size_t &len, const size_t &stride){
		double mean = 0.0;

		size_t nSteps = len/stride; // floor is what I want

		for (size_t i = 0; i < nSteps; i++) {
			mean += (arr[i * stride] - mean)/static_cast<double>(i + 1);
		}
		return mean;

	}

	/** \brief Mean of a C++ vector
	 *
	 * Calculates a mean of a vector.  Uses the recursive algorithm for numerical stability.
	 *
	 * \param[in] vec vector to average
	 *
	 */
	double mean(const vector<double> &vec){
		double mean = 0.0;

		for (size_t i = 0; i < vec.size(); i++) {
			mean += (vec[i] - mean)/static_cast<double>(i + 1);
		}
		return mean;
	}
	/** \brief Mean of a C++ vector with stride
	 *
	 * Calculates a mean of a vector with stride (i.e., using every _stride_ element).  Uses the recursive algorithm for numerical stability.
	 *
	 * \param[in] vec vector to average
	 * \param[in] stride stride
	 *
	 */
	double mean(const vector<double> &vec, const size_t &stride){
		double mean = 0.0;

		size_t nSteps = vec.size()/stride; // floor is what I want

		for (size_t i = 0; i < nSteps; i++) {
			mean += (vec[i * stride] - mean)/static_cast<double>(i + 1);
		}
		return mean;
	}

	/** \brief Square of a double
	 *
	 * \param[in] x value to square
	 * \return double square of the input
	 */
	inline double pow2(const double &x) {return x*x; };

	/** \brief Shift three values left
	 *
	 * Shifts a new value, moving the old to the left. The first value is discarded.
	 *
	 * \param[in,out] a first value (becomes _b_)
	 * \param[in,out] b second value (becomes _c_)
	 * \param[in,out] c third value (becomes _d_)
	 * \param[in] d unchanged fourth value (shifted to _c_)
	 *
	 */
	inline void shft3(double &a, double &b, double &c, const double &d){
		a = b;
		b = c;
		c = d;
	}

	/** \brief Bracket a maximum
	 *
	 * Brackets a maximum of a function given two initial guesses. Based on the Numerical Recipes in C++ function. Using max rather than min because max is more important in statistical applications.
	 * The resulting bracketing values are candA < candB < candC.
	 *
	 * \param[in] func function to maximize. Has to take a double and return a double
	 * \param[in] startA starting value A
	 * \param[in] startB starting value B
	 * \param[out] candA first bracketing value
	 * \param[out] candB second bracketing value
	 * \param[out] candC third breacketing value
	 *
	 */
	template<class T>
	void bracketMax(T &func, const double &startA, const double &startB, double &candA, double &candB, double &candC);

	// body of the function has to be in the header file because of template stuff
	template<class T>
	void bracketMax(T &func, const double &startA, const double &startB, double &candA, double &candB, double &candC){
		// set up constants
		const double GOLD   = 1.618034; // default successive magnification ratio
		const double GLIMIT = 100.0;    // maximum magnification allowed for the parabolic fit
		const double TINY   = 1.0e-20;  // tiny value to assure no division by zero

		// set the initial values
		candA = startA;
		candB = startB;
		double fa = func(candA);
		double fb = func(candB);

		// test if the function changes value at all (or at least a little over epsilon)
		if (fabs(fa - fb) < 1.001*BS_EPS) {
			candB += 100.0; // maybe just unlucky; change candB
			fb = func(candB);
			if (fabs(fa - fb) < 1.001*BS_EPS) { // now for sure there is a problem
				throw string("ERROR: function does not change over the candidate interval");
			}
		}
		// want to keep going uphill, so reverse order if b takes us down
		if (fb < fa) {
			swap(candA, candB);
			swap(fa, fb);
		}
		// first guess at candC
		candC = candB + GOLD*(candB - candA);

		double fc = func(candC);

		while (fb < fc) { // stop when we bracket (func(candC) drops below func(candB)); will not execute if our first guess at candC hit the jackpot

			// parabolic extrapolation from a,b,c to find a new candidate u
			// note that parabolic extrapolation works the same for min and max
			double r = (candB-candA)*(fb-fc);
			double q = (candB-candC)*(fb-fa);
			// NOTE: signbit() is C++11; true if negative
			double qrDiff = (signbit(q-r) ? -fmax(fabs(q-r), TINY) : fmax(fabs(q-r), TINY)); // using TINY to guard against division by zero
			double u      = candB - ((candB-candC)*q-(candB-candA)*r)/(2.0*qrDiff);
			double ulim   = candB + GLIMIT*(candC-candB);
			double fu;

			if (!signbit((candB - u)*(u - candC))) { // u is between b and c; try it
				fu = func(u);
				if (fu > fc) { // maximum between b and c
					candA = candB;
					candB = u;
					break;
				} else if (fu < fb){ // maximum between a and u
					candC = u;
					break;
				}
				// Nothing good found; revert to golden rule magnification
				u  = candC + GOLD*(candC-candB);
				fu = func(u);
			} else if (!signbit((candC-u)*(u-ulim))){ // u is between c and limit
				fu = func(u);
				if (fu > fc) {
					shft3(candB, candC, u, u+GOLD*(u-candC));
					shft3(fb, fc, fu, func(u)); // u has changed in the previous shift
				}
			} else if ((u-ulim)*(ulim-candC) >= 0.0) { // limit u to maximum allowed value if it is closer than c
				u  = ulim;
				fu = func(u);
			} else { // reject the parabolic approach and use golden rule magnification
				u  = candC + GOLD*(candC-candB);
				fu = func(u);
			}
			// did not bracket yet; eliminate the oldest point and continue
			shft3(candA, candB, candC, u);
			shft3(fa, fb, fc, fu);
		}

	}

	/** \brief Find the value that maximizes a function
	 *
	 * Uses the Brent method to find the value of \f$x\f$ that maximizes a function. Modification of the implementation found in Numerical Recipes in C++. Maximizing rather than minimizing because that is the most common application in statistics.
	 * Tolerance is set at \f$1.001 \times \sqrt{\epsilon}\f$, where \f$\epsilon\f$ is machine floating-point precision for _double_. This is just above the theoretical limit of precision.
	 *
	 * \param[in] func functor that represents the function to be maximized
	 * \param[in] startX starting value
	 * \param[out] xMax value of \f$x\f$ at maximum
	 * \param[out] fMax function value at maximum
	 *
	 */
	template<class T>
	void maximizer(T &func, const double &startX, double &xMax, double &fMax);


	// Template function. Body has to be in the .hpp file.
	template<class T>
	void maximizer(T &func, const double &startX, double &xMax, double &fMax){
		const double tol         = 1.001 * sqrt(BS_EPS); // set the tolerance just above the theoretical limit
		const unsigned int ITMAX = 1000;                 // maximum number of iterations
		const double CGOLD       = 0.3819660;            // golden ratio for when we abandon parabolic interpolation
		const double ZEPS        = 1e-3 * BS_EPS;        // small number to protect against numerical problems when the maximum is zero and we are trying to achieve a certain fractional accuracy

		// misc. parameters; named the same as the ones in Numerical Recipes Chapter 10.3
		double d = 0.0;
		double e = 0.0; // the distance moved in the step before last
		double etemp;
		double fu;
		double fv;
		double fw;
		double fx;
		double u;
		double v;
		double w;
		double x;
		double xm;
		double p;
		double q;
		double r;
		double tol1;
		double mtol1;
		double tol2;

		// start by bracketing
		const double startB = startX + 100.0;
		double ax;
		double bx;
		double cx;

		bracketMax(func, startX, startB, ax, bx, cx);

		double a = (ax < cx ? ax : cx);
		double b = (ax > cx ? ax : cx); // a and b (but not necessarily func(a) and func(b)) must be ascending order

		// initialize
		x  = w = v = bx;
		fw = fv = fx = func(x);
		unsigned int iter;
		for (iter = 0; iter < ITMAX; iter++) {
			xm    = 0.5 * (a + b); // xm is the midpoint between a and b
			tol1  = tol * fabs(x) + ZEPS;
			tol2  = 2.0 * tol1;
			mtol1 = -tol1;

			// test for doneness
			if (fabs(x - xm) <= (tol2 - 0.5*(b-a))) {
				fMax = fx;
				xMax = x;
				break;
			}

			if (fabs(e) > tol1) { // construct parabolic fit; not happening for the first round since e is set to 0.0 to begin with
				r = (x-w)*(fx-fv);
				q = (x-v)*(fx-fw);
				p = (x-v)*q-(x-w)*r;
				q = 2.0*(q-r);
				if (q > 0.0) {
					p = -p;
				}
				q     = fabs(q);
				etemp = e;
				e     = d;

				// Test the parabolic fit for acceptability.
				// The parabolic step has to be (1) in the (a,b) interval and (2) the movement has to be smaller than 0.5*(the one before last)
				if ((fabs(p) >= fabs(0.5*q*etemp)) || (p <= q*(a-x)) || (p >= q*(b-x))) {
					// parabolic step no good; take golden section into the larger segment
					e = (x >= xm ? a-x : b-x);
					d = CGOLD*e;
				} else { // take the parabolic step
					d = p/q;
					u = x+d;
					if ((u-a < tol2) || (b-u < tol2)) {
						d = (signbit(xm - x) ? mtol1 : tol1);
					}
				}
			} else {
				e = (x >= xm ? a-x : b-x);
				d = CGOLD * e;
			}
			u  = (fabs(d) >= tol1 ? x+d : x+(signbit(d) ? mtol1 : tol1));
			fu = func(u); // our one function evaluation

			// once we have our function evaluation, we decide what to do with it
			if (fu >= fx) {
				if (u >= x) {
					a = x;
				} else {
					b = x;
				}
				shft3(v,w,x,u);
				shft3(fv,fw,fx,fu);
			} else {
				if (u < x) {
					a = u;
				} else {
					b = u;
				}
				if ((fu >= fw) || (w == x)) {
					v  = w;
					w  = u;
					fv = fw;
					fw = fu;
				} else if ((fu >= fv) || (v == x) || (v == w)) {
					v  = u;
					fv = fu;
				}
			}
		}
		// if we did not get there after max # of iterations
		if (iter + 1 >= ITMAX) {
			xMax = x;
			fMax = nan("");

		}
	}

	/** \brief Logarithm of the Gamma function
	 *
	 * The log of the \f$ \Gamma(x) \f$ function. Implementing the Lanczos algorythm following Numerical Recipes in C++.
	 *
	 * \param[in] x value
	 * \return \f$ \log \Gamma(x) \f$
	 *
	 */
	double lnGamma(const double &x){
		if (x <= 0.0) return nan("");

		// define the weird magical coefficients
		const double coeff[14] {57.1562356658629235,-59.5979603554754912,14.1360979747417471,-0.491913816097620199,0.339946499848118887e-4,0.465236289270485756e-4,-0.983744753048795646e-4,0.158088703224912494e-3,-0.210264441724104883e-3,0.217439618115212643e-3,-0.164318106536763890e-3,0.844182239838527433e-4,-0.261908384015814087e-4,0.368991826595316234e-5};
		// save a copy of x for incrementing
		double y     = x;
		double gamma = 5.24218750000000000; // 671/128
		double tmp   = x + gamma;
		tmp          = (x + 0.5)*log(tmp) - tmp;
		double logPi = 0.91893853320467267;  // 0.5*log(2.0*pi)
		tmp         += logPi;
		double cZero = 0.999999999999997092; // c_0

		for (size_t i = 0; i < 14; i++) {
			cZero += coeff[i]/(++y);
		}

		return tmp + log(cZero/x);
	}

	/** \brief Continued fraction of the Beta function
	 *
	 * Computes the continued fraction of the Beta function following the Lenz method (see Numerical Recipes in C++). To be used in the _betai_ function.
	 *
	 * \param[in] x value
	 * \param[in] a shape parameter \f$a\f$
	 * \param[in] b shape parameter \f$b\f$
	 *
	 * \return continued fraction value
	 *
	 */
	double betacf(const double &x, const double &a, const double &b){
		const unsigned int maxIter = 10000;

		const double aPb = a + b;
		const double aPo = a + 1.0;
		const double aMo = a - 1.0;

		double numer = 1.0;             // first step of Lentz's method; define first continued fraction numerator
		double denom = 1.0 - aPb*x/aPo; // define the denominator to start

		if (fabs(denom) < BS_FPMIN) denom = BS_FPMIN; // set d to something resonalbly small if it is too close to 0
		denom = 1.0/denom;
		double cFrac = denom; // will become the continued fraction

		for (unsigned int m = 1; m < maxIter; m++) {
			double dm = static_cast<double>(m);
			double m2 = 2.0*dm;
			double aa = dm*(b-dm)*x/( (aMo+m2)*(a+m2) );

			// even step of the recurrence
			denom = 1.0 + aa*denom;
			numer = 1.0 + aa/numer;
			if (fabs(denom) < BS_FPMIN) denom = BS_FPMIN;
			if (fabs(numer) < BS_FPMIN) numer = BS_FPMIN;

			denom  = 1.0/denom;
			cFrac *= denom*numer;
			aa = -(a+dm)*(aPb+dm)*x/( (a+m2)*(aPo+m2) );

			// odd step of the recurrence
			denom = 1.0 + aa*denom;
			numer = 1.0 + aa/numer;
			if (fabs(denom) < BS_FPMIN) denom = BS_FPMIN;
			if (fabs(numer) < BS_FPMIN) numer = BS_FPMIN;

			denom = 1.0/denom;
			double del = denom*numer;
			cFrac *= del;
			if (fabs(del-1.0) <= BS_EPS) break; // done if within epsilon

		}
		return cFrac;
	}

	/** \brief Regularized incomplete Beta function
	 *
	 * Computes a quadrature approximatino of the regularized incomplete Beta function following the method in Numerical Recipes in C++. To be used in the _betai_ function.
	 *
	 * \param[in] x value
	 * \param[in] a shape parameter \f$a\f$
	 * \param[in] b shape parameter \f$b\f$
	 *
	 * \return approximate \f$ I_x(a, b) \f$
	 *s
	 */
	double betaiapprox(const double &x, const double &a, const double &b){
		// Gauss-Legendre abscissas and weights. Magic numbers copied from Numerical Recipes in C++
		const double y[18] {0.0021695375159141994,0.011413521097787704,0.027972308950302116,0.051727015600492421,0.082502225484340941, 0.12007019910960293,
			0.16415283300752470,0.21442376986779355, 0.27051082840644336,0.33199876341447887,0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
			0.62232745288031077,0.70331500465597174, 0.78649910768313447,0.87126389619061517, 0.95698180152629142};
		const double w[18] {0.0055657196642445571,0.012915947284065419,0.020181515297735382,0.027298621498568734,0.034213810770299537,0.040875750923643261,
			0.047235083490265582,0.053244713977759692,0.058860144245324798,0.064039797355015485,0.068745323835736408,0.072941885005653087,0.076598410645870640,
			0.079687828912071670,0.082187266704339706,0.084078218979661945,0.085346685739338721,0.085983275670394821};
		double res;
		double xu;
		const double aMo   = a - 1.0;
		const double bMo   = b - 1.0;
		const double mu    = a/(a+b);
		const double lnMu  = log(mu);
		const double lnOMU = log(1.0 - mu);

		double t = sqrt(a*b/(pow2(a+b)*(a+b+1.0)));

		// set the extent of tail integration
		if (x > mu) {
			if (x > 1.0) return 1.0;
			xu = fmax(mu + 10.0*t, x + 5.0*t);
			xu = fmin(1.0, xu);
		} else {
			if (x < 0.0) return 0.0;
			xu = fmin(mu - 10.0*t, x - 5.0*t);
			xu = fmax(0.0, xu);
		}
		double sum = 0.0;
		// Gauss-Legendre accumulation
		for (size_t i = 0; i < 18; i++) {
			t    = x + (xu-x)*y[i];
			sum += w[i]*exp(aMo*(log(t)-lnMu)+bMo*(log(1-t)-lnOMU));
		}
		res = sum*(xu-x)*exp(aMo*lnMu-lnGamma(a)+bMo*lnOMU-lnGamma(b)+lnGamma(a+b));
		return res > 0.0 ? 1.0-res : -res;
	}

	/** \brief Regularized incomplete Beta function
	 *
	 * Computes the regularized incomplete Beta function following the method in Numerical Recipes in C++.
	 *
	 * \param[in] x value
	 * \param[in] a shape parameter \f$a\f$
	 * \param[in] b shape parameter \f$b\f$
	 *
	 * \return \f$ I_x(a, b) \f$
	 *
	 */
	double betai(const double &x, const double &a, const double &b){
		if ( (a <= 0.0) || (b <= 0.0) ) return nan("");
		if ( (x < 0.0)  || (x > 1.0) )  return nan("");
		if ( (x == 0.0) || (x == 1.0) ) return x;

		const double doApprox = 3000.0; // when to do the quadrature approximation
		if ( (a > doApprox) && (b > doApprox) )  return betaiapprox(x, a, b);
		double bt = exp(lnGamma(a+b) - lnGamma(a) - lnGamma(b) + a*log(x) + b*log(1.0 - x));
		if (x < (a+1.0)/(a+b+2.0)) {
			return bt*betacf(x,a,b)/a;
		} else {
			return 1.0 - bt*betacf(1.0-x,b,a)/b;
		}
	}

	/** \brief Shell sort
	 *
	 * Sorts the provided vector in ascending order using Shell's method. Rather than move the elements themselves, save their indexes to the output vector. The first element of the index vector points to the smallest element of the input vector etc. The implementation is modified from code in Numerical Recipes in C++.
	 * NOTE: This algorithm is too slow for vectors of \f$ > 50\f$ elements. I am using it to finish off the quickSort, this is why I am giving it a range within a larger vector.
	 *
	 * \param[in] target vector to be sorted
	 * \param[in] beg index of the first element
	 * \param[in] end index of one past the last element to be included
	 * \param[out] outIdx vector of indexes
	 */
	void shellSort(const vector<double> &target, const size_t &beg, const size_t &end, vector<size_t> &outIdx){
		if (target.size() < end) {
			throw string("Target vector size smaller than end index in shellSort()");
		} else if (outIdx.size() < end) {
			throw string("Output vector size smaller than end index in shellSort()");
		} else if (end < beg) {
			throw string("End index smaller than beginning index in shellSort()");
		} else if (target.size() != outIdx.size()) {
			throw string("Target and output vectors must be of the same size in shellSort()");
		}
		// set up the initial output index values
		//for (size_t i = beg; i < end; i++) {
		//	outIdx[i] = i;
		//}
		// pick the initial increment
		size_t inc = 1;
		do {
			inc = inc*3 + 1;
		} while (inc <= end - beg);

		// start the sort
		do { // loop over partial sorts, decreasing the increment each time
			inc /= 3;
			const size_t bottom = beg + inc;
			for (size_t iOuter = bottom; iOuter < end; iOuter++) { // outer loop of the insertion sort, going over the indexes
				if (outIdx[iOuter] >= target.size()) {
					throw string("outIdx value out of bounds for target vector in shellSort()");
				}
				const size_t curInd = outIdx[iOuter]; // save the current value of the index
				size_t jInner       = iOuter;
				while (target[ outIdx[jInner - inc] ] > target[ curInd ]) {  // Straight insertion inner loop; looking for a place to insert the current value
					if (outIdx[jInner-inc] >= target.size()) {
						throw string("outIdx value out of bounds for target vector in shellSort()");
					}
					outIdx[jInner] = outIdx[jInner-inc];
					jInner        -= inc;
					if (jInner < bottom) {
						break;
					}
				}
				outIdx[jInner] = curInd;
			}
		} while (inc > 1);
	}
	/** Quicksort
	 *
	 * This function implements the Quicksort algorithm, taking the Numerical Recipes implementation as a base. It re-arranges the indexes in the output vector rather than move around the elements of the target vector. The output index must be the same size as the target (this is checked and exception thown if the condidition is not met). The output index is initialized with the correct index values.
	 *
	 * \param[in] target vector to be sorted
	 * \param[in] beg index of the first element
	 * \param[in] end index of one past the last element to be included
	 * \param[out] outIdx vector of indexes
	 *
	 */
	void quickSort(const vector<double> &target, const size_t &beg, const size_t &end, vector<size_t> &outIdx){
		if (target.size() < end) {
			throw string("Target vector size smaller than end index in quickSort()");
		} else if (outIdx.size() < end) {
			throw string("Output vector size smaller than end index in quickSort()");
		} else if (end < beg) {
			throw string("End index smaller than beginning index in quickSort()");
		} else if (target.size() != outIdx.size()) {
			throw string("Target and output vectors must be of the same size in quickSort()");
		}

		for (size_t i = beg; i < end; i++) {
			outIdx[i] = i;
		}
		// stacks to keep the l and ir values of the sub-vectors that are not being worked on
		stack<size_t> lStack;
		stack<size_t> irStack;

		const size_t m = 15; // the size of the sub-vectors that will be sorted by shellSort

		size_t l  = beg;     // left side of the sub-vector to be bisected
		size_t ir = end - 1; // right side of the sub-vector to be bisected

		while(true){
			if (ir-l < m) { // the current sub-vector is small enough for shellSort
				shellSort(target, l, ir+1, outIdx);
				if (lStack.empty()) {
					break;
				}
				l = lStack.top();
				lStack.pop();
				ir = irStack.top();
				irStack.pop();
			} else {
				const size_t k = (l+ir) >> 1; // median between l and ir
				swapXOR(outIdx[k], outIdx[l+1]);
				// rearrange the vector region so that target[outIdx[l]] <= target[outIdx[l+1]] <= target[outIdx[ir]]
				if (target[ outIdx[l] ] > target[ outIdx[ir] ]) {
					swapXOR(outIdx[l], outIdx[ir]);
				}
				if (target[ outIdx[l+1] ] > target[ outIdx[ir] ]) {
					swapXOR(outIdx[l+1], outIdx[ir]);
				}
				if (target[ outIdx[l] ] > target[ outIdx[l+1] ]) {
					swapXOR(outIdx[l], outIdx[l+1]);
				}

				// set the range for partitioning
				size_t i  = l+1;
				size_t j  = ir;
				size_t ip = outIdx[l+1]; // index of the pivot (partitioning element)
				while(true){ // inner loop
					do {
						i++;
					} while(target[ outIdx[i] ] < target[ ip ]); // scan forward to find element > pivot

					do {
						j--;
					} while(target[ outIdx[j] ] > target[ ip ]); // scan backwards to find element < pivot
					if (j < i) { // stop the scans if the indexes crossed
						break;
					}
					swapXOR(outIdx[i], outIdx[j]); // exchange elements unless the indexes have crossed
				}

				outIdx[l+1] = outIdx[j]; // insert the index of the pivot
				outIdx[j]   = ip;

				// push the indexes of defining the larger sub-vector to the stacks; the smaller sub-vector will be processed in the next iteration of the loop
				if ( (ir-i+1) > (j-l) ) {
					lStack.push(i);
					irStack.push(ir);
					ir = j - 1;
				} else {
					lStack.push(l);
					irStack.push(j-1);
					l = i;
				}
			}
		}
	}
}
#endif /* utilities_hpp */



