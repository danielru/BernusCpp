#include "bernus.h"
#include <vector>

// default constructor: initializes all gates to their steady-state value for V = -90.272 mV
bernus::bernus(std::vector<double>* gates):Iionmodel(gates) {

  (*this->gates).resize(bernus::ngates);

  double const Vrest = -90.272;
  
  (*this->gates)[m_gate]  = bnf.alpha_m(Vrest)/( bnf.alpha_m(Vrest) + bnf.beta_m(Vrest) );
  (*this->gates)[v_gate]  = bnf.v_inf(Vrest);
  (*this->gates)[f_gate]  = bnf.alpha_f(Vrest)/( bnf.alpha_f(Vrest) + bnf.beta_f(Vrest) );
  (*this->gates)[to_gate] = bnf.alpha_to(Vrest)/( bnf.alpha_to(Vrest) + bnf.beta_to(Vrest) );
  (*this->gates)[x_gate]  = bnf.x_inf(Vrest);
  
  gates_dt.resize(bernus::ngates);
  for (int i=0; i<bernus::ngates; ++i) {
    gates_dt[i] = 0.0;
  }
}

// destructor
bernus::~bernus() {
  // nothing to do
}