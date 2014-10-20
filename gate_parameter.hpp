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

#ifndef GATE_PARAMETER
#define GATE_PARAMETER

#include <cmath>

class gate_parameter
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
 * (1) Sodium current i_Na
 */

//! m-gate
inline double gate_parameter::alpha_m(double V)
{ return 0.32*(V+47.13)/(1.0 - exp(-0.1*(V+47.13))); }

inline double gate_parameter::beta_m(double V)
{ return 0.08*exp(-V/11.0); }


//! v-gate
inline double gate_parameter::v_inf(double V)
{ return 0.5*(1.0 - (tanh(7.74 + 0.12*V))); }

inline double gate_parameter::tau_v(double V)
{ return 0.25 + 2.24*( 1.0-(tanh(7.74 + 0.12*V)) )/( 1.0 - tanh(0.07*(92.4+V)) ); }

/**
 * (2) Calcium current i_Ca
 */

//! d-gate
inline double gate_parameter::d_inf(double V)
{ return alpha_d(V)/(alpha_d(V)+beta_d(V)); }

inline double gate_parameter::alpha_d(double V)
{ return 14.98*exp(-0.5)*pow(((V-22.36)/16.68), 2.0)/(16.68*sqrt(2.0*M_PI)); }

inline double gate_parameter::beta_d(double V)
{ return 0.1471 - (5.3*exp(-0.5)*pow( (V-6.27)/14.93, 2.0 ) )/(14.93*sqrt(2.0*M_PI)) ; }

//! f-gate
inline double gate_parameter::alpha_f(double V)
{ return 6.87*1e-3/(1.0 + exp( -(6.1546-V)/6.12) ); }

inline double gate_parameter::beta_f(double V)
{ return 5.75*1e-4 + (0.069*exp(-11.0*(V+9.825))+0.011)/(1.0 + exp(-0.278*(V+9.825))); }

//! f_Ca-gate
inline double gate_parameter::f_ca(double V)
{ return 1.0/(1.0 + this->Ca_i/0.0006); }

/**
 * (3) Transient outward current i_to
 */

//! r-gate
inline double gate_parameter::r_inf(double V)
{ return alpha_r(V)/(alpha_r(V)+beta_r(V)); }

inline double gate_parameter::alpha_r(double V)
{ return 0.5266*exp(-0.0166*(V-42.2912))/(1.0 + exp(-0.0943*(V-42.2912))); } //TODO: Insert correct function

inline double gate_parameter::beta_r(double V)
{ return (5.186*1e-5*V+0.5149*exp(-0.1344*(V-5.0027)))/(1.0 + exp(-0.1348*(V-5.186*1e-5))); } //TODO: Insert correct function

//! to-gate
inline double gate_parameter::alpha_to(double V)
{ return 0.0; } //TODO: Insert correct function

inline double gate_parameter::beta_to(double V)
{ return 0.0; } //TODO: Insert correct function

inline double gate_parameter::tau_to(double V)
{ return 0.0; } //TODO: Insert correct function

inline double gate_parameter::to_inf(double V)
{ return 0.0; } //TODO: Insert correct function

/**
 * (4) Delayed rectifier potassium current i_K
 */

//! X-gate
inline double gate_parameter::x_inf(double V)
{ return 0.0; } //TODO: Insert correct function

inline double gate_parameter::tau_x(double V)
{ return 0.0; } //TODO: Insert correct function

inline double gate_parameter::tau_x_a(double V)
{ return 0.0; } //TODO: Insert correct function

/**
 * (5) Inward rectifier potassium current i_K1
 */

//! K1-gate
inline double k1_inf(double V)
{ return 0.0; } //TODO: Insert correct function

inline double alpha_k1(double V)
{ return 0.0; } //TODO: Insert correct function

inline double beta_k1(double V)
{ return 0.0; } //TODO: Insert correct function

/**
 * (8) Sodium potassium pump
 */
inline double f_nak(double V)
{
double sigma = 0.0;
return 0.0; } //TODO: Insert correct function

inline double f_nak_a(double V)
{ return 0.0; } //TODO: Insert correct function

/**
 * (9) Sodium calcium pump i_NaCa
 */
inline double f_naca(double V)
{ return 0.0; } //TODO: Insert correct function

#endif // GATE_PARAMETER