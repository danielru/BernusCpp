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
  
  std::cout << "m_t:   " << brn.gates_dt[brn.m_gate] << std::endl;
  std::cout << "v_t:   " << brn.gates_dt[brn.v_gate] << std::endl;
  std::cout << "f_t:   " << brn.gates_dt[brn.f_gate] << std::endl;
  std::cout << "to_t:  " << brn.gates_dt[brn.to_gate] << std::endl;
  std::cout << "x_t:   " << brn.gates_dt[brn.x_gate] << std::endl;
  std::cout << "Iion:  " << Iion << std::endl;

}