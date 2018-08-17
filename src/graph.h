#pragma once

#include <vector>

class Graph
{
	std::vector<double> unary_costs;
	std::vector<double> pairwise_costs;

public:

	int node_count, edge_count, stride;

	Graph();
	~Graph();

	void load_from_file(std::string path);
	
};