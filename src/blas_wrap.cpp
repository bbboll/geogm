#include "blas_wrap.h"
#include <assert.h>
#include <iostream>
#include <cstring>

#ifdef __APPLE__
	#include <Accelerate/Accelerate.h>
#else
	#include <cblas.h>
#endif

/*
	Numerical comparison.
*/
bool almost_equal(double a, double b) {
	return std::abs( a - b ) < EPS;
}

/*
	pointwise subdivide by operand
	B ← A / B
*/
void matrix_psubdivide(const Matrix& A, Matrix& B) {
	assert( A.m == B.m );
	assert( A.n == B.n );

	// the following will usually compile down 
	// to vector instructions on x86_64
	double *dataA = A.raw();
	double *dataB = B.raw();
	for (int i = 0; i < A.m*A.n; ++i)
	{
		dataB[i] = dataA[i] / dataB[i];
	}
}

/*
	pointwise exponential
	A ← exp(A)
*/
void matrix_pexp(Matrix &A) {
	double *dataA = A.raw();
	for (int i = 0; i < A.m*A.n; ++i)
	{
		dataA[i] = std::exp(dataA[i]);
	}
}

/*
	pointwise natural logarithm
	A ← log(A)
*/
void matrix_plog(Matrix &A) {
	double *dataA = A.raw();
	for (int i = 0; i < A.m*A.n; ++i)
	{
		dataA[i] = std::log(dataA[i]);
	}
}

/*
	pointwise power
	A ← A**p
*/
void matrix_ppower(Matrix &A, const double p) {
	double *dataA = A.raw();
	for (int i = 0; i < A.m*A.n; ++i)
	{
		dataA[i] = std::pow(dataA[i], p);
	}
}

/*
	pointwise product
	A ← A*B
*/
void matrix_pproduct(Matrix &A, const Matrix& B) {
	double *dataA = A.raw();
	double *dataB = B.raw();
	for (int i = 0; i < A.m*A.n; ++i)
	{
		dataA[i] = dataA[i]*dataB[i];
	}
}

/*
	Mean squared error
*/
double matrix_mse(Matrix &A, Matrix &B) {
	assert( A.m == B.m );
	assert( A.n == B.n );

	double *dataA = A.raw();
	double *dataB = B.raw();
	double sum = 0.0;
	for (int i = 0; i < A.m*A.n; ++i)
	{
		sum += (dataA[i] - dataB[i]) * (dataA[i] - dataB[i]);
	}
	return sum / (A.m*A.n);
}

/*
	BLAS routine daxpy
	x ← x + αy
*/
void matrix_add(Matrix &x, const double alpha, const Matrix&y) {
	assert( x.m*x.n == y.m*y.n );

	cblas_daxpy(
		x.m*x.n,	// number of entries
		alpha,      // multiplier alpha
		y.raw(),	// other data
		1,			// other stride
		x.raw(),	// x data
		1			// x stride
	);
}

/*
	BLAS routine copy
	B ← A
*/
void matrix_copy(const Matrix& A, Matrix& B) {
	assert( A.m == B.m );
	assert( A.n == B.n );

	cblas_dcopy(
		A.m*A.n,	// number of entries
		A.raw(),	// source data
		1,			// source stride
		B.raw(),	// destination data
		1			// destination stride
	);
}

/*
	BLAS routine dgemv
	y ← αAx + βy
*/
void matrix_vector_multiply(const Matrix &A, const Matrix &x, Matrix& y, double alpha, double beta) {
	auto transA = CBLAS_TRANSPOSE::CblasNoTrans;
	if (A.trans) {
		transA = CBLAS_TRANSPOSE::CblasTrans;
	}

	assert( A.second_dim() == x.m );
	assert( A.first_dim() == y.m );
	assert( x.n == 1 );
	assert( y.n == 1 );

	cblas_dgemv(
		CBLAS_ORDER::CblasRowMajor, // memory layout
		transA,                     // transpose A
		A.m,                        // left   dimension
		A.n,                        // right  dimension
		alpha,                      // alpha multiplier
		A.raw(),                    // A data
		A.n,                        // lda
		x.raw(),                    // x data
		1,                          // x stride
		beta,                       // beta multiplier
		y.raw(),                    // y data
		1                           // y stride
	);
}

/*
	BLAS routine dgemm
	C ← αAB + βC
*/
void matrix_multiply(const Matrix& A, const Matrix& B, Matrix& C, double alpha, double beta) {

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

	// DGEMM routine from accelerate framework or openBLAS (cblas interface)
	cblas_dgemm(
		CBLAS_ORDER::CblasRowMajor, // memory layout
		transA,                     // transpose A
		transB,                     // transpose B
		C.m,                        // left   dimension
		C.n,                        // right  dimension
		A.second_dim(),             // middle dimension
		alpha,                      // alpha multiplier
		A.raw(),                    // A data
		A.n,                        // lda
		B.raw(),                    // B data
		B.n,                        // ldb
		beta,                       // beta multiplier
		C.raw(),                    // C data
		C.n                         // ldc
	);
}

/*
	BLAS routine ddot
*/
double dot_product(const Matrix& A, const Matrix& B) {
	assert( A.m*A.n == B.m*B.n );

	return cblas_ddot(
		A.m*A.n,   // N dimension
		A.raw(),   // A data
		1,         // A stride
		B.raw(),   // B data
		1          // B stride
	);
}

Matrix& Matrix::operator+=(const Matrix& other) {
	matrix_add(*this, 1.0, other);
	return *this;
};