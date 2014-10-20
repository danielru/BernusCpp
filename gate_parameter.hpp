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
  double alpha_d(double);
  double beta_d(double);
  double alpha_f(double);
  double beta_f(double);
  double f_ca(double);
  
  //! Transient outward current
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
  double sigma(double);
  
  //! Sodium calcium pump
  double f_naca(double);
  
};

/**
 * Implementation of class functions; kept in header for easier inlining
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
{ return 0.0; } //TODO: Insert correct function

inline double gate_parameter::tau_v(double V)
{ return 0.0; } //TODO: Insert correct function

/**
 * (2) Calcium current i_Ca
 */

//! d-gate
inline double gate_parameter::alpha_d(double V)
{ return 0.0; } //TODO: Insert correct function

inline double gate_parameter::beta_d(double V)
{ return 0.0; } //TODO: Insert correct function

//! f-gate
inline double gate_parameter::alpha_f(double V)
{ return 0.0; } //TODO: Insert correct function

inline double gate_parameter::beta_f(double V)
{ return 0.0; } //TODO: Insert correct function

//! f_Ca-gate
inline double gate_parameter::f_ca(double V)
{ return 0.0; } //TODO: Insert correct function

/**
 * (3) Transient outward current i_to
 */

//! r-gate
inline double gate_parameter::alpha_r(double V)
{ return 0.0; } //TODO: Insert correct function

inline double gate_parameter::beta_r(double V)
{ return 0.0; } //TODO: Insert correct function

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

/**
 * (6) Calcium background current i_b_Ca
 */

/**
 * (7) Sodium background current i_b_Na
 */

/**
 * (8) Sodium potassium pump
 */

/**
 * (9) Sodium calcium pump i_NaCa
 */

#endif // GATE_PARAMETER