#include "graph.h"
#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

void Graph::load_from_file(std::string path) {
	// Get decompressed stream from file
	std::ifstream file(path, std::ios_base::in | std::ios_base::binary);
	boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
	inbuf.push(boost::iostreams::gzip_decompressor());
	inbuf.push(file);
	std::istream instream(&inbuf);

	std::string line;
	while(std::getline(instream, line)) {
	    std::cout << line << std::endl;
	}
	
	file.close();
}