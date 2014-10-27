#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <time.h>
#include "Iionmodel.h"
#include "IionmodelFactory.h"

int main(int args, char** argv) {
  
  std::fstream output_file;

  double const capacitance = 1.0;
  double V0   = -60;
  double Tend = 500.0;
  int nsteps  = 1e7;
  double dt   = Tend/( (double) nsteps );
  double Iion;
  clock_t timer = clock();
  bool output = false;
  
  std::cout << "Time step (ms): " << dt << std::endl;
  
  std::vector<double> gates;
  std::vector<double> gates_dt;
  
  Iionmodel * brn = IionmodelFactory::factory(IionmodelFactory::bernus, &gates, &gates_dt);
  
  output_file.open("./bernus.txt", std::ios_base::out);
  
  for(int i=0; i<nsteps; ++i) {
    
    if ( (i<100) || (i % 1000 == 0) ) {
      output = true;
      output_file << dt*( (double) i) << "    ";
      output_file << V0 << "    ";
    }
    else{ output = false;}
    
    // Update derivative of gating variables
    brn->update_gates_dt(V0);
    
    // Compute ionic currents
    Iion = brn->ionforcing(V0);
    
    // Forward Euler update for gating variables
    for (int j=0; j<brn->get_ngates(); ++j) {
      (*brn->gates)[j] += dt*(*brn->gates_dt)[j];
      
      if (output) {
        output_file << (*brn->gates)[j] << "    ";
      }
    }
    
    // Forward Euler update of membrane potential
    V0 += -(1.0/capacitance)*dt*Iion;
    
    if (output) {
      output_file << Iion << std::endl;
    }
  }
  
  output_file.close();
  
  timer = clock() - timer;
  float time_in_sec = ( (float) timer )/CLOCKS_PER_SEC;
  std::cout << "Total runtime:                       " << time_in_sec << std::endl;
  std::cout << "Average time per ion model timestep: " << time_in_sec/( (double) nsteps) << std::endl;

  // Print out steady-state values for gating variables:
  std::cout << std::endl;
  return 0;
}