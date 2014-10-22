/**
 * NOTE: Multiple small one-line functions are implemented in the header file bernus.hpp
 */

#include "bernus.hpp"
#include <vector>

// constructor
bernus::bernus():Iionmodel() {
  gates.resize(bernus::ngates);
  for(size_t i=0; i<bernus::ngates; ++i){
    gates[i] = 1.0;
  }
  gates_dt.resize(bernus::ngates);
}

// destructor
bernus::~bernus() {
  // bla
}
