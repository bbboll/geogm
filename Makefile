CXX=clang++
SRC=src/*.cpp
STD=c++1z

blastest: src/blastest.cpp src/blas_wrap.h
	$(CXX) -framework Accelerate \
		   -std=$(STD) \
		   -o blastest src/blastest.cpp

read_gz: src/read_gz.cpp
	$(CXX) -std=$(STD) \
		   -lboost_iostreams \
		   -o read_gz src/read_gz.cpp

clean:
	rm blastest