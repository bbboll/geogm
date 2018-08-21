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

/*
	Check for zero entries in a matrix and
	rectify them.
*/
void rectify(Matrix& A) {
	const double EPS = 1e-10;
	double min = 1.0;
	double* dataA = A.raw();
	for (int i = 0; i < A.m*A.n; ++i)
	{
		if (dataA[i] < min)
		{
			min = dataA[i];
		}
	}
	if (min < EPS)
	{
		A += (double) (- min + EPS);
		A *= (double) (1 / A.sum());
	}
}

// wasserstein gradient iteration limit
// and convergence threshold
const int ITERATION_LIMIT = 15;
const double CONVERGENCE_THRESH = 1e-6;
const int FIRST_COMPONENT = 1;
const int SECOND_COMPONENT = 2;

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
		const double tau,
		const int component,
		Matrix& gradient_component) {

	// stride
	const int n = Theta.n;
	 
	// allocate iteration variables
	// TODO: eliminate allocations by using a single memory pool for all threads
	Matrix v_current(n, 1);
	Matrix K(n, n);
	Matrix _tmp(n, 1);

	// init iteration variables
	v_current.one_fill();

	// init K
	matrix_copy(Theta, K);
	K *= (-1/tau);
	matrix_pexp(K); 

	for (int i = 0; i < ITERATION_LIMIT; ++i)
	{
		if (component == 1) {
			// v1 update
			K.transpose();
			matrix_vector_multiply(K, v_current, _tmp, 1.0, 0.0);
			matrix_psubdivide(mu2, _tmp);
			K.transpose();
			matrix_vector_multiply(K, _tmp, gradient_component, 1.0, 0.0);
			matrix_psubdivide(mu1, gradient_component);
		}
		else {
			// v2 update
			matrix_vector_multiply(K, v_current, _tmp, 1.0, 0.0);
			matrix_psubdivide(mu2, _tmp);
			K.transpose();
			matrix_vector_multiply(K, _tmp, gradient_component, 1.0, 0.0);
			matrix_psubdivide(mu2, gradient_component);
		}

		if ( matrix_mse(v_current, gradient_component) < CONVERGENCE_THRESH ) {
			break;
		}
	}

	// recover optimal dual values
	matrix_plog(gradient_component);
	gradient_component *= tau;

	// project onto tangent space
	project_T(gradient_component);
}

/*
	4.13
*/
void energy_gradient_row(
		const Graph& g, 
		const int i, 
		const double tau, 
		Matrix& energy_gradient) {
	std::vector<int> outbound = g.adjacent_edges(i).second;
	std::vector<int> inbound = g.adjacent_edges(i).second;

	matrix_copy(g.get_unary(i), energy_gradient);
	project_T(energy_gradient);

	Matrix component(g.stride, 1);

	for (auto const& ij: outbound) {
		wasserstein_gradient(
			g.get_unary_labels(g.edge_endpoints.at(ij).first), 
			g.get_unary_labels(g.edge_endpoints.at(ij).second), 
			g.get_pairwise(ij), 
			tau, 
			FIRST_COMPONENT, 
			component
		);
		energy_gradient += component;
	}
	for (auto const& ji: inbound) {
		wasserstein_gradient(
			g.get_unary_labels(g.edge_endpoints.at(ji).first), 
			g.get_unary_labels(g.edge_endpoints.at(ji).second), 
			g.get_pairwise(ji), 
			tau, 
			SECOND_COMPONENT, 
			component
		);
		energy_gradient += component;
	}
}

/*
	6.1
*/
void iterate(
		const Graph& g, 
		const double tau, 
		const double alpha, 
		const double h) {
	// setup iteration variable memory
	Matrix energy_gradient(g.stride, 1);
	Matrix unary_labels(g.stride, 1);

	for (int i = 0; i < g.node_count; ++i)
	{
		energy_gradient_row(g, i, tau, energy_gradient);
		energy_gradient *= (-h);
		matrix_pexp(energy_gradient);

		matrix_copy(g.get_unary_labels(i), unary_labels);
		matrix_ppower(unary_labels, 1 + alpha);

		const double denom = dot_product(unary_labels, energy_gradient);
		matrix_pproduct(unary_labels, energy_gradient);
		unary_labels *= (1/denom);
		rectify(unary_labels);
		Matrix update_labels = g.get_unary_labels(i);
		matrix_copy(unary_labels, update_labels);
	}

}

int main( int argc, char **argv )
{
	if(argc < 2) {
		std::cerr << "Usage: optimize <gzipped input file>" << std::endl;
		return 0;
	}

	Graph g;
	g.load_from_file(argv[1]);

	// setup hyperparameters
	const double tau   = 0.1;
	const double alpha = 1.5;
	const double h     = 0.1;
	const int GLOBAL_ITERATION_LIMIT = 15;

	for (int i = 0; i < GLOBAL_ITERATION_LIMIT; ++i)
	{
		std::cout << g.energy() << std::endl;
		iterate(g, tau, alpha, h);
	}

	std::cout << "Final energy: " << g.energy() << std::endl;

	return 0;
}