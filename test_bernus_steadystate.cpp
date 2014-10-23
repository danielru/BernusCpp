/**
 * Verifies that Iion(Vrest) = 0 for Vrest = -90.272 mV
 */
#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <time.h>
#include "bernus.hpp"
#include "bernus.cpp"

int main(int args, char** argv) {
  
  double const Vrest = -90.272;
  bernus brn;

  //brn.plot_equil_potentials();
  
  brn.update_gates_dt(Vrest);
  double Iion = brn.ionforcing(Vrest);
  
  for (int i=0; i<brn.ngates; ++i){
    std::cout << brn.gates_dt[i] << std::endl;
  }
  
  std::cout << "Iion: " << Iion << std::endl;

}