CXX=clang++
SRC=src/*.cpp
STD=c++1z

blastest: src/blastest.cpp src/blas_wrap.h
	$(CXX) -framework Accelerate \
		   -std=$(STD) \
		   -o blastest src/blastest.cpp

clean:
	rm blastest