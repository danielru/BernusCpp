#include "bernus.hpp"
#include "bernus.cpp"
#include <cstdlib>
#include <iostream>
#include <fstream>

int main(int args, char** argv)
{
  
  std::fstream output_file;
  
  double V = -90, V_max = 50;
  int    steps = 1400;
  double V_inc = std::abs(V_max - V)/( (double) steps );
  
  double mgate;
  double vgate;
  double fgate;
  double togate;
  double xgate;
  double Iion;
  
  bernus brn;
  
  output_file.open("./bernus_probe.txt", std::ios_base::out);
  
  for (int i=0; i<steps; ++i) {
    
    
    //V += V_inc;
    //output_file << V << "    " << v_inf << "    " << tau_v << "    " << x_inf << "     " << tau_x << std::endl;
  }
  
  output_file.close();
  return 0;
}