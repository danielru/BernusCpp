#include "bernus_functions.hpp"

int main(int args, char** argv)
{
  double V = 0.01;
  int ncalls = 1e8;
  bernus_functions * bf = new bernus_functions();

  double alpha_m, beta_m, v_inf, sum;
  sum = 0.0;
  for(int i=0; i<ncalls; i++)
  {
    alpha_m = bf->alpha_m(V);
    beta_m  = bf->beta_m(V);
    v_inf   = bf->v_inf(V);
    sum = sum + alpha_m + beta_m + v_inf;

  }
  delete bf;
  return 0;
}