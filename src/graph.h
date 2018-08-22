#pragma once

#include "blas_wrap.h"
#include <vector>
#include <utility>
#include <string>

typedef std::pair<std::vector<int>,std::vector<int>> intvecpair;

class Graph
{
	double *unary_costs = NULL;
	double *pairwise_costs = NULL;
	double *unary_labels = NULL;
	std::vector< intvecpair > adjacency;
	Matrix unary_costs_m;
	Matrix pairwise_costs_m;
	Matrix unary_labels_m;

public:

	int node_count, edge_count, stride;
	std::vector< std::pair<int, int> > edge_endpoints;

	Graph() {};
	~Graph() {
		if (unary_costs != NULL) {
			free(unary_costs);
		}
		if (pairwise_costs != NULL) {
			free(pairwise_costs);
		}
		if (unary_labels != NULL) {
			free(unary_labels);
		}
	};

	void load_from_file(std::string path);
	const Matrix get_unary(const int i) const;
	const Matrix get_pairwise(const int ij) const;
	Matrix get_unary_labels(const int i) const;
	const intvecpair adjacent_edges(const int i) const;
	const double energy() const;
	void export_labeling(const int width, const int height) const;
	
};