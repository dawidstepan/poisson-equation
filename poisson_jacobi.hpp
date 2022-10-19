#pragma once

#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#include <omp.h>
#include <assert.h>
#include <Eigen/Dense>

using namespace Eigen;

template <typename Type> class MatrixView {
private:
  std::vector<Type> &v;
  MatrixView(const MatrixView &);
  MatrixView &operator=(const MatrixView &);

public:
  const size_t N, M;
  MatrixView(std::vector<Type> &v, size_t N, size_t M) : v(v), N(N), M(M) {
    assert(v.size() / N == M);
  }
  Type &set(size_t i, size_t j) { return v[i + N * j]; }
  const Type &get(size_t i, size_t j) { return v[i + N * j]; }
  Type &set(size_t n) { return v[n]; }
  const Type &get(size_t n) { return v[n]; }
};

double ParticularSolution(double x, double y) {
  return sin(2 * M_PI * x) * sinh(2 * M_PI * y);
}

double NormL2(const std::vector<double> &v) {
  double norm = 0;
  for (const auto &value : v) {
    norm += value * value;
  }
  return sqrt(norm);
}

double NormInf(const std::vector<double> &v) {
  double max = std::numeric_limits<double>::lowest();
  for (const auto &value : v) {
    max = std::fabs(value) > max ? std::fabs(value) : max;
  }
  return max;
}

struct Stencil {
  Stencil(double h)
      : C(4.0 / (h * h) + 4 * M_PI * M_PI), N(-1.0 / (h * h)),
        S(-1.0 / (h * h)), W(-1.0 / (h * h)), E(-1.0 / (h * h)) {}
  const double C, N, S, W, E;
};

enum Cell { UNKNOWN = 0, DIR = 1, NEU = 2, ROB = 0 };


void PoissonJacobiStencil(size_t resolution,size_t threads,int *iterations,double *runtime, double *runtime2) {
 
  
  size_t NY = resolution;
  size_t NX = (2 * NY)-1;
  double h = 1.0 / (NY - 1);


  const auto stencil = Stencil(h);

  // domain cell types
  std::vector<int> domain(NX * NY, Cell::UNKNOWN);
  MatrixView<int> domainView(domain, NX, NY);
  for (size_t i = 0; i != NX; ++i) {
    domainView.set(i, 0) = Cell::DIR;
    domainView.set(i, NY - 1) = Cell::DIR;
  }
  for (size_t j = 0; j != NY; ++j) {
    domainView.set(0, j) = Cell::DIR;
    domainView.set(NX - 1, j) = Cell::DIR;
  }

  // referenceSolution
  std::vector<double> referenceSolution(NX * NY, 0);
  MatrixView<double> referenceSolutionView(referenceSolution, NX, NY);
  for (size_t j = 0; j < NY; ++j) {
    for (size_t i = 0; i < NX; ++i) {
      referenceSolutionView.set(i, j) = ParticularSolution(i * h, j * h);
    }
  }

  // right hand side
  std::vector<double> rightHandSide(NX * NY, 0);
  MatrixView<double> rightHandSideView(rightHandSide, NX, NY);
  for (size_t j = 0; j < NY; ++j) {
    for (size_t i = 0; i < NX; ++i) {
      rightHandSideView.set(i, j) =
          ParticularSolution(i * h, j * h) * 4 * M_PI * M_PI;
    }
  }

  auto SolverJacobi = [&](std::vector<double> &sol, std::vector<double> &sol2,
                         std::vector<double> &rhs, const Stencil &stencil,
                         size_t NX, size_t NY) {
    MatrixView<double> solView(sol, NX, NY);
    MatrixView<double> sol2View(sol2, NX, NY);
    MatrixView<double> rhsView(rhs, NX, NY);

   // Set threads and parallize
   omp_set_num_threads(threads); 
   #pragma omp parallel for collapse(2) 

    for (size_t j = 1; j < NY - 1; ++j) {
       for (size_t i = 1; i < NX - 1; ++i) {

        sol2View.set(i, j) =
            1.0 / stencil.C *
            (rhsView.set(i, j) - (solView.get(i + 1, j) * stencil.E +
                                  solView.get(i - 1, j) * stencil.W +
                                  solView.get(i, j + 1) * stencil.S +
                                  solView.get(i, j - 1) * stencil.N));
      }
    }
    sol.swap(sol2);
  };

  // solution approximation starting with boundary initialized to dirichlet
  // conditions, else 0
  std::vector<double> solution(NX * NY, 0);
  MatrixView<double> solutionView(solution, NX, NY);
  for (size_t j = 0; j != NY; ++j) {
    for (size_t i = 0; i != NX; ++i) {
      if (domainView.get(i, j) == Cell::DIR)
        solutionView.set(i, j) = ParticularSolution(i * h, j * h);
    }
  }
  std::vector<double> solution2 = solution;
  std::cout << "Solve LSE using stencil jacobi" << std::endl;
  auto start = std::chrono::high_resolution_clock::now();
  
  int secondsin{10};
  int num{0}; 
  double runtimein{0};

  while (std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::system_clock::now() - start)
	     .count() < secondsin) {
			  ++num;
  		          auto startin = std::chrono::high_resolution_clock::now();
			  SolverJacobi(solution, solution2, rightHandSide, stencil, NX, NY);
			  auto stopin = std::chrono::high_resolution_clock::now();
			  auto secondsin =
      			  std::chrono::duration_cast<std::chrono::duration<double>>(stopin - startin)
          		      .count();
			  runtimein += secondsin;			  

  }

  auto stop = std::chrono::high_resolution_clock::now();
  auto seconds =
       std::chrono::duration_cast<std::chrono::duration<double>>(stop - start)
           .count();
  

  std::cout << "Iterations=" << num << std::endl;
  std::cout << "Runtime=" << std::scientific << seconds << std::endl;
  std::cout << "Average runtime per iteration=" << seconds/num << std::endl;
  std::cout << "Average runtime per iteration (within)=" << runtimein/num << std::endl;
	
  *iterations=num;
  *runtime=seconds/num;
  *runtime2=runtimein/num;
 
}
