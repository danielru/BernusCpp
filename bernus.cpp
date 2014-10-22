/**
 * NOTE: Multiple small one-line functions are implemented in the header file bernus.hpp
 */

#include "bernus.hpp"

// constructor
bernus::bernus():Iionmodel() {
  gates.resize(this->ngates);
  
  //TODO: Somehow initialize gating variables
  for(size_t i=0; i<this->ngates; ++i){
    gates[i] = 0.0;
  }
}

// destructor
bernus::~bernus() {
  // bla
}
