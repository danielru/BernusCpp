#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <time.h>
#include "Iionmodel.h"
#include "IionmodelFactory.h"

// For testing
#include "bernus.h"

int main(int args, char** argv) {
  
  std::fstream output_file;

  double const capacitance = 1.0;
  double V0   = -20;
  double Tend = 3;
  int nsteps  = 500;
  double dt   = Tend/( (double) nsteps );
  double Iion;
  clock_t timer = clock();
  bool output = false;
  bool repol = false;
  
  std::cout << "Time step (ms): " << dt << std::endl;
  
  std::vector<double> gates;
  std::vector<double> gates_dt;
  
  Iionmodel * brn = IionmodelFactory::factory(IionmodelFactory::BERNUS, &gates, &gates_dt);
  
  // For testing, want to also access functions in Bernus which are not exposed by the interface
  bernus * bbb = (bernus*) brn;
  
  output_file.open("./bernus.txt", std::ios_base::out);
  
  for(int i=0; i<nsteps; ++i) {
    
    if ( (i<250) || (i % 1 == 0) ) {
      output = true;
      output_file << dt*( (double) i) << "    ";
      output_file << V0 << "    ";
    }
    else{ output = false;}
    
    // Compute ionic currents
    Iion = brn->ionforcing(V0);
    
    // Rush-Larsen update of gates
    brn->rush_larsen_step(V0, dt);
    
    if (output) {
      for (int j=0; j<brn->get_ngates(); ++j) {
        output_file << (*brn->gates)[j] << "    ";
      }
    }
    
    if ( (i*dt>25.0) && !repol && (brn->gates[bernus::m_gate]<0.98)) {
      std::cout << "Repolarized at t = " << i*dt << std::endl;
      repol = true;
    }
    
    // Forward Euler update of membrane potential
    V0 += -(1.0/capacitance)*dt*Iion;
    
    if (output) {
      //output_file << Iion << std::endl;
      
      output_file << (bbb->i_na)(V0) << std::endl;
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