#include "bernus_functions.hpp"
#include <cstdlib>
#include <iostream>
#include <fstream>

int main(int args, char** argv)
{
  
  std::fstream output_file;
  
  double V = -90, V_max = 50;
  int    steps = 1400;
  double V_inc = std::abs(V_max - V)/( (double) steps );
  
  double v_inf, tau_v, x_inf, tau_x;
  
  bernus_functions brn;
  
  output_file.open("./bernus_functions.txt", std::ios_base::out);
  
  for (int i=0; i<steps; ++i) {
    v_inf = brn.v_inf(V);
    tau_v = brn.tau_v(V);
    x_inf = brn.x_inf(V);
    tau_x = brn.tau_x(V);
    V += V_inc;
    output_file << V << "    " << v_inf << "    " << tau_v << "    " << x_inf << "     " << tau_x << std::endl;
  }
  
  output_file.close();
  return 0;
}