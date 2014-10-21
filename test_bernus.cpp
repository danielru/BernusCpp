#include <cstdlib>
#include "bernus.hpp"
#include <omp.h>
#include <iostream>

int main(int narg, char** args)
{
  int nruns = 1e8;
  bernus * bnm = new bernus();
  double timer = omp_get_wtime();
  double dummy = 0.0;
  
  for(int i=0; i<nruns; i++)
  {
    double a = (double) rand() / RAND_MAX;
    dummy += bnm->ionforcing(a);
  }

  timer = omp_get_wtime() - timer;

  std::cout << "Total runtime:                    " << timer << std::endl;
  std::cout << "Average time per ionforcing call: " << timer/( (double) nruns) << std::endl;
}