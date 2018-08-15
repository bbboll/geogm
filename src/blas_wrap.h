#ifndef BLAS_WRAPPER_HH
#define BLAS_WRAPPER_HH

#include <assert.h>
#include <iostream>
#include <cmath>

#ifdef __APPLE__
	#include <Accelerate/Accelerate.h>
#else
	#include <cblas.h>
#endif

const double EPS = 1e-14;

/*
	Numerical comparison.
*/
bool almost_equal(double a, double b) {
	return abs( a - b ) < EPS;
}

/*
	Type to be used with BLAS wrappers.
	This also hides C-style memory management.
*/
class Matrix
{
	double *A;

public:
	// matrix dimensions
	const int m, n;

	// transposed?
	bool trans = false;

	Matrix(const int _m, const int _n) : m(_m), n(_n) {
		// allocate memory for entries
		A = (double*) malloc( m*n * sizeof(double) );
	};

	~Matrix() { 
		free(A);
	};

	double *raw() {
		return A;
	};

	// toggle transpose state
	void transpose() {
		trans = !trans;
	};

	int second_dim() const {
		if (trans) {
			return m;
		}
		return n;
	};

	int first_dim() const {
		if (trans) {
			return n;
		}
		return m;
	};

	double get(int i, int j) const {
		assert( i < first_dim() );
		assert( j < second_dim() );

		if (trans) {
			return A[ j*n + i ];
		}
		return A[ i*n + j ];
	};

	void set(int i, int j, double x) {
		assert( i < first_dim() );
		assert( j < second_dim() );

		if (trans) {
			A[ j*n + i ] = x;
		}
		else {
			A[ i*n + j ] = x;
		}
	};

	void zero_fill() {
		memset( A, 0.0, m*n*sizeof(double) );
	};

	void one_fill() {
		for (int i = 0; i < first_dim(); ++i)
		{
			for (int j = 0; j < second_dim(); ++j)
			{
				set(i,j, 1.0);
			}
		}
	};

	void print() const {
		std::cout << "Matrix " << first_dim() << "x" << second_dim() << std::endl;
		for (int i = 0; i < first_dim(); ++i)
		{
			for (int j = 0; j < second_dim(); ++j)
			{
				std::cout << " " << get(i,j);
			}
			std::cout << std::endl;
		}
	};
};

/*
	BLAS routine dgemm
	C ← αAB + βC
*/
void matrix_multiply(Matrix& A, Matrix& B, Matrix& C, double alpha, double beta) {

	auto transA = CBLAS_TRANSPOSE::CblasNoTrans;
	if (A.trans) {
		transA = CBLAS_TRANSPOSE::CblasTrans;
	}
	auto transB = CBLAS_TRANSPOSE::CblasNoTrans;
	if (B.trans) {
		transB = CBLAS_TRANSPOSE::CblasTrans;
	}

	assert( A.second_dim() == B.first_dim() );
	assert( A.first_dim()  == C.m );
	assert( B.second_dim() == C.n );

	// DGEMM routine from accelerate framework (cblas interface)
	cblas_dgemm(
		CBLAS_ORDER::CblasRowMajor, // memory layout
		transA,                     // transpose A
		transB,                     // transpose B
		C.m,                        // left   dimension
		C.n,                        // right  dimension
		A.second_dim(),             // middle dimension
		alpha,                      // alpha multiplier
		A.raw(),                    // A data
		A.n,                        // lda >= left_dim
		B.raw(),                    // B data
		B.n,                        // ldb >= right_dim
		beta,                       // beta multiplier
		C.raw(),                    // C data
		C.n                         // ldc >= right_dim
	);
}

#endif