#include <assert.h>
#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>
#include "arguments.hpp"
#include "poisson_jacobi.hpp"

 
int main(int argc, char *argv[]) {

  // Parse command line arguments
  auto resolution = convertTo<int>(1, 32, argc, argv);
  auto threads = convertTo<int>(2, 800, argc, argv);
  auto schedule = argv[3]; //for printing issues

  assert(resolution > 0);
  assert(threads > 0);

  std::cout << "Input:" << std::endl;
  std::cout << "resolution=" << resolution << std::endl;
  std::cout << "threads=" << threads << std::endl;
  std::cout << "schedule=" << schedule << std::endl;

  int iterations{0};
  double runtime{0};
  double runtime2{0};

  std::cout << "Output:" << std::endl;
  PoissonJacobiStencil(resolution,threads,&iterations,&runtime,&runtime2);

  // Test output 
  std::cout << "Max. number of available threads: " << omp_get_max_threads() << std::endl;

  // Print results to txt-file
  std::ofstream log;
  log.open("results.txt",std::ofstream::app);
  log << resolution << "," << threads << "," << schedule << "," << runtime << "," << runtime2 << std::endl;
  cout << std::endl;

  return 0;

}
