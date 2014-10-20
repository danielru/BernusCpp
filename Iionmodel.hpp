#ifndef IIONMODEL
#define IIONMODEL

class Iionmodel {

public:

  Iionmodel() {};
  
  // destructor declared virtual to ensure proper polymorphic delete
  virtual ~Iionmodel() {};
  
  double ionforcing(double);
  
  double* statevars_rhs();
  
};


#endif // IIONMODEL