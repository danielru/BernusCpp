#ifndef IIONMODEL
#define IIONMODEL

#include <vector>

class Iionmodel {

public:

  Iionmodel() {};
  
  // destructor declared virtual to ensure proper polymorphic delete
  virtual ~Iionmodel() {};
  
  virtual double ionforcing(double) = 0;
  
  virtual std::vector<double> statevars_rhs(double) = 0;
  
};


#endif // IIONMODEL