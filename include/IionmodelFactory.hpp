#include "Iionmodel.hpp"
#include "bernus.hpp"
#include "bernus.cpp"
#include <stdexcept>

/**
 * @todo docu
 */
class IionmodelFactory {
  
public:
  
  //! @todo docu
  static const int bernus = 1;
  
  //! @todo docu
  static Iionmodel * factory(int choice) {
    if (choice == IionmodelFactory::bernus)
      return bernus::factory();
    else
      throw std::runtime_error("IIonmodelFactory: No model available for selected value of choice");
  }
  
};