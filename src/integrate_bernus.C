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
  double const Vrest = -90.272;
  double const Vmax = 0;
  double const eps_recovery    = 0.01;
  double const eps_development = 0.04;
  double const kTa  = 1;
  double const gamma = 0.5;           // Newmark integrator parameter

  double const Vpace = 30.272;
  double V0;
  double const Tend = 500;
  int const nsteps  = 1e4;
  int const npacing = 2;
  double const dt   = Tend/( (double) nsteps );
  double Iion;
  double Ta = 0;
  clock_t timer = clock();
  bool output = false;
  bool repol = false;
  
  std::cout << "Time step (ms): " << dt << std::endl;
  
  std::vector<double> gates;
  
  Iionmodel * brn = IionmodelFactory::factory(IionmodelFactory::BERNUS);
  
  brn->initialize(&gates);
  
  
  // For testing, want to also access functions in Bernus which are not exposed by the interface
  bernus * bbb = (bernus*) brn;
  
  output_file.open("./bernus.txt", std::ios_base::out);
  
  for(int npace=0; npace<npacing; npace++) {
      V0 = Vrest + Vpace;
      
	  for(int i=0; i<nsteps; ++i) {
	
		if ( (i<250) || (i % 1 == 0) ) {
		  output = true;
		  output_file << dt*( (double) i)+npace*Tend << "    ";
		  output_file << V0 << "    ";
		}
		else{ output = false;}
	
		// Compute ionic currents
		Iion = brn->ionforcing(V0, &gates);
	
		// Rush-Larsen update of gates
		brn->rush_larsen_step(V0, dt, &gates);
	
		if (output) {
		  for (int j=0; j<brn->get_ngates(); ++j) {
			output_file << gates[j] << "    ";
		  }
		}
	
		if ( (i*dt>25.0) && !repol && (gates[bernus::m_gate]<0.98)) {
		  std::cout << "Repolarized at t = " << i*dt << std::endl;
		  std::cout << "Potential at this time = " << V0 << std::endl;
		  repol = true;
		}
	
	    double const Vo = V0;
		// Forward Euler update of membrane potential
		V0 += -(1.0/capacitance)*dt*Iion;
	
		// normalized potential
		double const V0n = (V0 - Vrest)/(Vmax-Vrest);
		double const Von = (Vo - Vrest)/(Vmax-Vrest);
		// shape function
		auto eps = [&] (double V) { return (V < 0.05 ? eps_development : eps_recovery); };
		// Extended mean value theorem / Newmark integrator
		Ta =    (dt*((1-gamma)*eps(Von)*(kTa*Von-Ta)+gamma*eps(V0n)*kTa*V0n)+Ta)
		      / (1+dt*gamma*eps(V0n));
	
		if (output) {
		  //output_file << Iion << std::endl;
	  
		  output_file << (bbb->i_na)(V0,&gates) << "   " << Ta << "   " << V0n << std::endl;
		}
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