UNAME_S := $(shell uname -s)
BLAS_DIR := /usr/local

CXX=clang++
SRC=src/*.cpp
STD=c++1z

ifeq ($(UNAME_S),Linux)
blastest: src/blastest.cpp src/blas_wrap.h
		$(CXX) -I$(BLAS_DIR)/include \
			   -lopenblas \
			   -lc \
			   -std=$(STD) \
			   -o blastest src/blastest.cpp
endif
ifeq ($(UNAME_S),Darwin)
blastest: src/blastest.cpp src/blas_wrap.h
	$(CXX) -framework Accelerate \
		   -std=$(STD) \
		   -o blastest src/blastest.cpp
endif

read_gz: src/read_gz.cpp
	$(CXX) -std=$(STD) \
		   -lboost_iostreams \
		   -o read_gz src/read_gz.cpp

optimize: src/optimize.cpp src/graph.cpp
	$(CXX) -std=$(STD) \
		   -lboost_iostreams \
		   -lboost_regex \
		   -o optimize src/optimize.cpp src/graph.cpp

clean:
	rm optimize
	rm blastest
	rm read_gz