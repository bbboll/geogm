#pragma once

#include <assert.h>
#include <iostream>
#include <cmath>
#include <cstring>

const double EPS = 1e-14;

/*
	Numerical comparison.
*/
bool almost_equal(double a, double b);

/*
	Type to be used with BLAS wrappers.
	This also hides underlying memory management.
*/
class Matrix
{
	double *A = NULL;
	bool needs_memory_cleanup = true;

public:
	// matrix dimensions
	int m, n;

	// transposed?
	bool trans = false;

	Matrix() : m(0), n(0) {
		needs_memory_cleanup = false;
	};

	Matrix(const int _m, const int _n) : m(_m), n(_n) {
		// allocate memory for entries
		A = (double*) malloc( m*n * sizeof(double) );
	};

	Matrix(const int _m, const int _n, double *_A) : m(_m), n(_n) {
		// use externally managed memory location for A
		A = _A;
		needs_memory_cleanup = false;
	};

	~Matrix() {
		if (needs_memory_cleanup && A != NULL) {
			free(A);
		}
	};

	void set_data(double* data) {
		assert( A == NULL );
		A = data;
	};

	void set_size(const int _m, const int _n) {
		assert( m == 0 );
		assert( n == 0 );
		m = _m;
		n = _n;
	};

	double* raw() const {
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

	Matrix& operator*=(const double alpha) {
		for (int i = 0; i < m*n; ++i)
		{
			A[i] *= alpha;
		}
		return *this;
	};

	Matrix& operator+=(const double alpha) {
		for (int i = 0; i < m*n; ++i)
		{
			A[i] += alpha;
		}
		return *this;
	};

	Matrix& operator+=(const Matrix& other);

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

	double sum() const {
		double sum = 0.0;
		for (int i = 0; i < m*n; ++i)
		{
			sum += A[i];
		}
		return sum;
	};

	int argmax() const {
		double max = 0.0;
		int arg = 0;
		for (int i = 0; i < m*n; ++i)
		{
			if (A[i] > max)
			{
				max = A[i];
				arg = i;
			}
		}
		return arg;
	}
};

/*
	pointwise subdivide by operand
	B ← A / B
*/
void matrix_psubdivide(const Matrix& A, Matrix& B);

/*
	pointwise exponential
	A ← exp(A)
*/
void matrix_pexp(Matrix &A);

/*
	pointwise natural logarithm
	A ← log(A)
*/
void matrix_plog(Matrix &A);

/*
	pointwise power
	A ← A**p
*/
void matrix_ppower(Matrix &A, const double p);

/*
	pointwise product
	A ← A*B
*/
void matrix_pproduct(Matrix &A, const Matrix& B);

/*
	Mean squared error
*/
double matrix_mse(Matrix &A, Matrix &B);

/*
	BLAS routine daxpy
	y ← y + αx
*/
void matrix_add(Matrix &y, const double alpha, const Matrix&x);

/*
	BLAS routine copy
	B ← A
*/
void matrix_copy(const Matrix& A, Matrix& B);

/*
	BLAS routine dgemv
	y ← αAx + βy
*/
void matrix_vector_multiply(const Matrix &A, const Matrix &x, Matrix& y, double alpha, double beta);

/*
	BLAS routine dgemm
	C ← αAB + βC
*/
void matrix_multiply(const Matrix& A, const Matrix& B, Matrix& C, double alpha, double beta);

/*
	BLAS routine ddot
*/
double dot_product(const Matrix& A, const Matrix& B);