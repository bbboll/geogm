#pragma once

#include <vector>
#include <string>

class Graph
{
	double *unary_costs = NULL;
	double *pairwise_costs = NULL;

public:

	int node_count, edge_count, stride;

	Graph() {};
	~Graph() {
		if (unary_costs != NULL) {
			free(unary_costs);
		}
		if (pairwise_costs != NULL) {
			free(pairwise_costs);
		}
	};

	void load_from_file(std::string path);
	
};