#include "blas_wrap.h"
#include "graph.h"
#include <assert.h>
#include <iostream>

/*
	Project a given vector onto the tangent space of the
	(n-1) dimensional probability simplex.
*/
void project_T(Matrix& x) {
	double mean = x.sum() / x.m;
	x += -mean;
}

// wasserstein gradient iteration limit
// and convergence threshold
const int ITERATION_LIMIT = 15;
const double CONVERGENCE_THRESH = 1e-6;

/*
	Local Wasserstein gradient.
	Maximize dual program of Wasserstein distance by Sinkhorn's algorithm. 
	Projecting the dual maximizers onto the tangent space yields the
	local Wasserstein gradient (Theorem).
*/
void wasserstein_gradient(
		const Matrix& mu1, 
		const Matrix& mu2, 
		const Matrix& Theta, 
		double tau) {

	// stride
	const int n = Theta.n;
	 
	// allocate iteration variables
	// TODO: eliminate allocations by using a single memory pool for all threads
	Matrix v1_current(n, 1);
	Matrix v1_next(n, 1);
	Matrix v2_current(n, 1);
	Matrix v2_next(n, 1);
	Matrix K(n, n);
	Matrix _tmp(n, 1);

	// init iteration variables
	v1_current.one_fill();
	v2_current.one_fill();

	// init K
	matrix_copy(Theta, K);
	K *= (-1/tau);
	matrix_pexp(K); 

	for (int i = 0; i < ITERATION_LIMIT; ++i)
	{
		// v1 update
		K.transpose();
		matrix_vector_multiply(K, v1_current, _tmp, 1.0, 0.0);
		matrix_psubdivide(mu2, _tmp);
		K.transpose();
		matrix_vector_multiply(K, _tmp, v1_next, 1.0, 0.0);
		matrix_psubdivide(mu1, v1_next);

		// v2 update
		matrix_vector_multiply(K, v2_current, _tmp, 1.0, 0.0);
		matrix_psubdivide(mu2, _tmp);
		K.transpose();
		matrix_vector_multiply(K, _tmp, v2_next, 1.0, 0.0);
		matrix_psubdivide(mu2, v2_next);

		// check convergence
		if (   matrix_mse(v1_current, v1_next) < CONVERGENCE_THRESH 
			&& matrix_mse(v2_current, v2_next) < CONVERGENCE_THRESH)
		{
			break;
		}
	}

	// recover optimal dual values
	matrix_plog(v1_next);
	matrix_plog(v2_next);
	v1_next *= tau;
	v2_next *= tau;

	// project onto tangent space
	project_T(v1_next);
	project_T(v2_next);
}

/*
	4.13
*/
void energy_gradient(const Graph& g, const int i) {
	std::vector<int> outbound = g.adjacent_edges(i).second;
	std::vector<int> inbound = g.adjacent_edges(i).second;

	for (auto const& ij: outbound) {
		// 
	}
	for (auto const& ij: inbound) {
		// 
	}
}

/*
	6.1
*/
void iterate(Matrix& W, Matrix& Theta, const double alpha, const double tau, const int n) {
	// TODO
}

int main( int argc, char **argv )
{
	if(argc < 2) {
		std::cerr << "Usage: optimize <gzipped input file>" << std::endl;
		return 0;
	}

	Graph g;
	g.load_from_file(argv[1]);

	// std::vector<int> v = g.adjacent_edges(1).second;
	// for (int i = 0; i < v.size(); ++i)
	// {
	// 	std::cout << v.at(i) << std::endl;	
	// }

	// energy_gradient(g, 1);
	// Matrix a = g.get_pairwise(1);
	// a.print();


	return 0;
}