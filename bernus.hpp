#ifndef BERNUS
#define BERNUS

#include "Iionmodel.hpp"
#include "bernus_functions.hpp"

class bernus: public Iionmodel {

public:
  
  bernus() {
    //TODO: Find out how to make this const and static but still a function of
    // the other parameters...
    this->test = bernus::R*bernus::T*bernus::Fa;
    for(int i=0; i<this->ngates; i++){
      this->gates[i] = 0.0;
    }
  };

  double ionforcing(double);
  
  double* statevars_rhs();
  
  //! The number of gating variables described by an ODE
  int static const ngates = 5;
  
private:

  
  /**
   * The total ionic current of the Bernus model is comprised of nine 
   * different currents.
   */
  
   //! Sodium current
   double i_na(double);
  
   //! Calcium current
   double i_ca(double);
  
   //! Transient outward current
   double i_to(double);
  
   //! Delayed rectifier potassium current
   double i_k(double);
  
  //! Inward rectifier potassium current
  double i_k1(double);
  
  //! Calcium background current
  double i_b_ca(double);
  
  //! Sodium background current
  double i_b_na(double);
  
  //! Sodium potassium pump
  double i_na_k(double);
  
  //! Sodium calcium pump
  double i_na_ca(double);
  
  /**
   * Member variables
   */
  
   //! Values of ODE-based gate variables
   double gates [ngates];

  /**
   * Constants below are from Table 1 in Bernus et al.
   * TODO: Check physical units...
   */
   //TODO: Should make precision here choosable as well
  double static const g_na   = 16.0;
  double static const g_ca   = 0.064;
  double static const g_to   = 0.4;
  double static const g_k    = 0.019;
  double static const g_k1   = 3.9;
  double static const g_na_b = 0.001;
  double static const g_ca_b = 0.00085;
  double static const g_nak  = 1.3;
  double static const g_naca = 1.0;
  
  //! Universal gas constant
  double static const R = 0.0; //TODO: Correct value
  
  //! Absolute temperature
  double static const T = 0.0; //TODO: Correct value
  
  //! Faraday constant
  double static const Fa = 0.0; //TODO: Correct value
  
  //! Extra- and inner cellular potentials (Table 1 in Bernus et al.)
  
  //! Equilibrium potentials
  double test;

};
#endif // BERNUS