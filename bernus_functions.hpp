/**
 * Header file that provides a collection of functions the
 * compute parameters of the  model for a given
 * membrane voltage potential.
 *
 * All functions here are oneliners, so that the compiler
 * can inline them easily.
 *
 * Daniel Ruprecht, October 20, 2014
 *
 */

#ifndef BERNUS_FUNCTIONS
#define BERNUS_FUNCTIONS

#include <cmath>

class bernus_functions
{
  
  public:
  
  //! TODO: To later be able to easily run the whole ion channel
  //! model in single precision, introduce a ionprec typedef accuracy
  
  //! Sodium current
  double alpha_m(double);
  double beta_m(double);
  double v_inf(double);
  double tau_v(double);
  
  //! Calcium current
  double d_inf(double);
  double alpha_d(double);
  double beta_d(double);
  double alpha_f(double);
  double beta_f(double);
  double f_ca(double);
  
  //! Transient outward current
  double r_inf(double);
  double alpha_r(double);
  double beta_r(double);
  double alpha_to(double);
  double beta_to(double);
  double tau_to(double);
  double to_inf(double);
  
  //! Delayed rectifier potassium current
  double x_inf(double);
  double tau_x(double);
  double tau_x_a(double);
  
  //! Inward rectifier potassium current
  double k1_inf(double);
  double alpha_k1(double);
  double beta_k1(double);
  
  //! Calcium background current: No parameter functions needed
  
  //! Sodium background current: No parameter functions needed
  
  //! Sodium potassium pump
  double f_nak(double);
  double f_nak_a(double);
  
  //! Sodium calcium pump
  double f_naca(double);
  
  private:
  
  double static const Ca_i = 0.0; // TODO: Insert correct value
  
};

/**
 * Implementation of class functions; kept in header for easier inlining.
 * See e.g. https://models.cellml.org/e/5/bernus_wilders_zemlin_verschelde_panfilov_2002.cellml/@@cellml_math
 * for the formulas.
 */

/**
 * (1) Sodium current i_Na (4 Functions)
 */

//! m-gate
inline double bernus_functions::alpha_m(double V)
{ return 0.32*(V+47.13)/(1.0 - exp(-0.1*(V+47.13))); }

inline double bernus_functions::beta_m(double V)
{ return 0.08*exp(-V/11.0); }

//! v-gate
inline double bernus_functions::v_inf(double V)
{ return 0.5*(1.0 - (tanh(7.74 + 0.12*V))); }

inline double bernus_functions::tau_v(double V)
{ return 0.25 + 2.24*( 1.0-(tanh(7.74 + 0.12*V)) )/( 1.0 - tanh(0.07*(92.4+V)) ); }

/**
 * (2) Calcium current i_Ca (5 functions)
 */

//! d-gate
inline double bernus_functions::d_inf(double V)
{ return alpha_d(V)/(alpha_d(V)+beta_d(V)); }

inline double bernus_functions::alpha_d(double V)
{ return 14.98*exp(-0.5)*pow(((V-22.36)/16.68), 2.0)/(16.68*sqrt(2.0*M_PI)); }

inline double bernus_functions::beta_d(double V)
{ return 0.1471 - (5.3*exp(-0.5)*pow( (V-6.27)/14.93, 2.0 ) )/(14.93*sqrt(2.0*M_PI)) ; }

//! f-gate
inline double bernus_functions::alpha_f(double V)
{ return 6.87*1e-3/(1.0 + exp( -(6.1546-V)/6.12) ); }

inline double bernus_functions::beta_f(double V)
{ return 5.75*1e-4 + (0.069*exp(-11.0*(V+9.825))+0.011)/(1.0 + exp(-0.278*(V+9.825))); }

//! f_Ca-gate
inline double bernus_functions::f_ca(double V)
{ return 1.0/(1.0 + this->Ca_i/0.0006); }

/**
 * (3) Transient outward current i_to (7 functions)
 */

//! r-gate
inline double bernus_functions::r_inf(double V)
{ return alpha_r(V)/(alpha_r(V)+beta_r(V)); }

inline double bernus_functions::alpha_r(double V)
{ return 0.5266*exp(-0.0166*(V-42.2912))/(1.0 + exp(-0.0943*(V-42.2912))); }

inline double bernus_functions::beta_r(double V)
{ return (5.186*1e-5*V+0.5149*exp(-0.1344*(V-5.0027)))/(1.0 + exp(-0.1348*(V-5.186*1e-5))); }

//! to-gate
inline double bernus_functions::alpha_to(double V)
{ return (5.612e-5*V+0.0721*exp(-0.173*(V+34.2531)))/(1.0 + exp(-0.1604*(V+34.0235))); } //TODO: Insert correct function

inline double bernus_functions::beta_to(double V)
{ return 0.0; } //TODO: Insert correct function

inline double bernus_functions::tau_to(double V)
{ return 0.0; } //TODO: Insert correct function

inline double bernus_functions::to_inf(double V)
{ return 0.0; } //TODO: Insert correct function

/**
 * (4) Delayed rectifier potassium current i_K (3 functions)
 */

//! X-gate
inline double bernus_functions::x_inf(double V)
{ return 0.0; } //TODO: Insert correct function

inline double bernus_functions::tau_x(double V)
{ return 0.0; } //TODO: Insert correct function

inline double bernus_functions::tau_x_a(double V)
{ return 0.0; } //TODO: Insert correct function

/**
 * (5) Inward rectifier potassium current i_K1 (3 functions)
 */

//! K1-gate
inline double k1_inf(double V)
{ return 0.0; } //TODO: Insert correct function

inline double alpha_k1(double V)
{ return 0.0; } //TODO: Insert correct function

inline double beta_k1(double V)
{ return 0.0; } //TODO: Insert correct function

/**
 * (8) Sodium potassium pump (3 functions)
 */
inline double f_nak(double V)
{
double sigma = 0.0;
return 0.0; } //TODO: Insert correct function

inline double f_nak_a(double V)
{ return 0.0; } //TODO: Insert correct function

/**
 * (9) Sodium calcium pump i_NaCa (1 function)
 */
inline double f_naca(double V)
{ return 0.0; } //TODO: Insert correct function

#endif // BERNUS_FUNCTIONS