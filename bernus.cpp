/**
 * NOTE: Multiple small one-line functions are implemented in the header file bernus.hpp
 */

#include "bernus.hpp"
#include <vector>

// default constructor: initializes all gates to their steady-state value for V = -90.272 mV
bernus::bernus():Iionmodel() {
  gates.resize(bernus::ngates);

  double const Vrest = -90.272;
  
  gates[m_gate]  = bnf.alpha_m(Vrest)/( bnf.alpha_m(Vrest) + bnf.beta_m(Vrest) );
  gates[v_gate]  = bnf.v_inf(Vrest);
  gates[f_gate]  = bnf.alpha_f(Vrest)/( bnf.alpha_f(Vrest) + bnf.beta_f(Vrest) );
  gates[to_gate] = bnf.alpha_to(Vrest)/( bnf.alpha_to(Vrest) + bnf.beta_to(Vrest) );
  gates[x_gate]  = bnf.x_inf(Vrest);
  
  gates_dt.resize(bernus::ngates);
  for (int i=0; i<bernus::ngates; ++i) {
    gates_dt[i] = 0.0;
  }
}

// destructor
bernus::~bernus() {
  // bla
}

void bernus::plot_equil_potentials(){
  std::cout << "Equilibrium potentials: " << std::endl;
  std::cout << "E_na: " << bnf.e_na << std::endl;
  std::cout << "E_ca: " << bnf.e_ca << std::endl;
  std::cout << "E_to: " << bnf.e_to  << std::endl;
  std::cout << "E_k: "  << bnf.e_k  << std::endl;
}