#include "graph.h"
#include <fstream>
#include <iostream>
#include <assert.h>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/regex.hpp>

void Graph::load_from_file(std::string path) {
	// build regex for parsing
	boost::regex all_ints{"(\\d+)"};
	boost::regex all_doubles{"(\\d+\\.\\d+)"};

	// Get decompressed stream from file
	std::ifstream file(path, std::ios_base::in | std::ios_base::binary);
	boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
	inbuf.push(boost::iostreams::gzip_decompressor());
	inbuf.push(file);
	std::istream instream(&inbuf);

	std::string line;

	// assert valid file header
	if (std::getline(instream, line)) {
		assert( line == "MARKOV" );
	}
	else {
		std::cout << "Encountered empty file." << std::endl;
		return;
	}

	// read number of nodes
	if (std::getline(instream, line)) {
		// read ints in line
		boost::regex_token_iterator<std::string::iterator> it{line.begin(), line.end(), all_ints, 1};
		boost::regex_token_iterator<std::string::iterator> end;
		while (it != end) {
			node_count = std::stoi(*it++);
		}
	}
	else {
		std::cout << "Illegal file format. Expected number of nodes but found EOF." << std::endl;
		return;
	}

	// read number of labels per node
	if (std::getline(instream, line)) {
		// read ints in line
		boost::regex_token_iterator<std::string::iterator> it{line.begin(), line.end(), all_ints, 1};
		boost::regex_token_iterator<std::string::iterator> end;
		int label_count = -1;
		while (it != end) {
			if (label_count == -1) {
				label_count = std::stoi(*it++);
			}
			else {
				// we only support constant label count
				assert( label_count == std::stoi(*it++) );
			}
		}
		stride = label_count;
	}
	else {
		std::cout << "Illegal file format. Expected number of labels per node but found EOF." << std::endl;
		return;
	}

	// read number of edges
	if (std::getline(instream, line)) {
		// read ints in line
		boost::regex_token_iterator<std::string::iterator> it{line.begin(), line.end(), all_ints, 1};
		boost::regex_token_iterator<std::string::iterator> end;
		while (it != end) {
			edge_count = std::stoi(*it++) - node_count;
		}
	}
	else {
		std::cout << "Illegal file format. Expected number of nodes but found EOF." << std::endl;
		return;
	}

	// skip node id assignment, we assume it to be ascending.
	for (int i = 0; i < node_count; ++i)
	{
		if (std::getline(instream, line)) {
			std::string prefix("1 ");
			auto mismatch_index = std::mismatch(prefix.begin(), prefix.end(), line.begin());
			assert( mismatch_index.first == prefix.end() );
		}
		else {
			std::cout << "Illegal file format. Expected node id assignment but found EOF." << std::endl;
		}
	}

	// read edge assignment
	for (int i = 0; i < edge_count; ++i)
	{
		if (!std::getline(instream, line)) {
			std::cout << "Illegal file format. Expected node-node assignment but found EOF." << std::endl;
		}
		// assert valid prefix
		std::string prefix("2 ");
		auto mismatch_index = std::mismatch(prefix.begin(), prefix.end(), line.begin());
		assert( mismatch_index.first == prefix.end() );

		// TODO -----------------------------------------------------------------------
	}

	// allocate memory for unary costs
	unary_costs = (double*) malloc( stride*node_count * sizeof(double) );
	int unary_cost_index = 0;

	// read unary costs
	for (int i = 0; i < node_count; ++i)
	{
		// first line is empty
		if (!std::getline(instream, line)) {
			std::cout << "Illegal file format. Expected unary cost but found EOF." << std::endl;
		}
		assert( line.size() == 0 );

		// second line contains stride
		if (!std::getline(instream, line)) {
			std::cout << "Illegal file format. Expected unary cost but found EOF." << std::endl;
		}
		boost::regex_token_iterator<std::string::iterator> ind_it{line.begin(), line.end(), all_ints, 1};
		boost::regex_token_iterator<std::string::iterator> end;
		while (ind_it != end) {
			assert( std::stoi(*ind_it++) == stride );
		}

		// third line contains cost
		if (!std::getline(instream, line)) {
			std::cout << "Illegal file format. Expected unary cost but found EOF." << std::endl;
		}
		boost::regex_token_iterator<std::string::iterator> double_it{line.begin(), line.end(), all_doubles, 1};
		while (double_it != end) {
			unary_costs[unary_cost_index++] = std::stod(*double_it++);
		}
	}

	assert( unary_cost_index == stride*node_count );

	// allocate memory for pairwise costs
	pairwise_costs = (double*) malloc( stride*stride*edge_count * sizeof(double) );
	int pairwise_cost_index = 0;

	// read pairwise costs
	for (int i = 0; i < edge_count; ++i)
	{
		// first line is empty
		if (!std::getline(instream, line)) {
			std::cout << "Illegal file format. Expected pairwise cost but found EOF." << std::endl;
		}
		assert( line.size() == 0 );

		// second line contains stride
		if (!std::getline(instream, line)) {
			std::cout << "Illegal file format. Expected pairwise cost but found EOF." << std::endl;
		}
		boost::regex_token_iterator<std::string::iterator> ind_it{line.begin(), line.end(), all_ints, 1};
		boost::regex_token_iterator<std::string::iterator> end;
		while (ind_it != end) {
			assert( std::stoi(*ind_it++) == stride*stride );
		}

		// third line contains cost
		if (!std::getline(instream, line)) {
			std::cout << "Illegal file format. Expected pairwise cost but found EOF." << std::endl;
		}
		boost::regex_token_iterator<std::string::iterator> double_it{line.begin(), line.end(), all_doubles, 1};
		while (double_it != end) {
			pairwise_costs[pairwise_cost_index++] = std::stod(*double_it++);
		}
	}

	assert( pairwise_cost_index == stride*stride*edge_count );
	
	file.close();
}