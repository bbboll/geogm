#include "graph.h"
#include <assert.h>
#include <iostream>

int main( int argc, char **argv )
{
	if(argc < 2) {
		std::cerr << "Usage: optimize <gzipped input file>" << std::endl;
		return 0;
	}

	Graph g;
	g.load_from_file(argv[1]);

	return 0;
}