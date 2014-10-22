#include <cstdlib>
#include "bernus.hpp"
#include "bernus.cpp"
#include <time.h>
#include <iostream>

int main(int narg, char** args)
{
  int nruns = 1e8;
  bernus * bnm = new bernus();
  clock_t timer = clock();
  double dummy = 0.0;
  
  for(int i=0; i<nruns; i++)
  {
    double a = (double) rand() / RAND_MAX;
    dummy += bnm->ionforcing(a);
  }

  timer = clock() - timer;

  float time_in_sec = ( (float) timer )/CLOCKS_PER_SEC;
  std::cout << "Total runtime:                    " << time_in_sec << std::endl;
  std::cout << "Average time per ionforcing call: " << time_in_sec/( (double) nruns) << std::endl;
}