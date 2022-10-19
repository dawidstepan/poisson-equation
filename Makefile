# requirements on ubuntu
# sudo apt-get build-essentials
CXX=g++

# Update the path to your Eigen library by replacing '/home/~'
EIGEN=/home/~

CXXFLAGS=-std=c++14 -O3 -Wall -pedantic -march=native -ffast-math -fopenmp

.DEFAULT_GOAL := all

all: stenciljacobi
 
stenciljacobi: Makefile solver.cpp poisson_jacobi.hpp arguments.hpp
	$(CXX) solver.cpp -o stenciljacobi -lpthread -v $(CXXFLAGS) -I$(EIGEN)

.PHONY: clean
clean:
	rm stenciljacobi
