#ifndef BERNUS_HPP
#define BERNUS_HPP
//Define NDEBUG to activate runtime asserts
#define NDEBUG
#include <algorithm>
#include <vector>
#include <assert.h>
#include "Iionmodel.h"
#include "bernus_functions.h"

/**
 * Class implementing the Bernus et al. model for ventricular cells:
 *
 * O. Bernus, R. Wilders, C. W. Zemlin, H. Verschelde, A. V. Panfilov:
 * "A computationally efficient electrophysiological model of human ventricular cells";
 * Am. J. Physiol. Heart Circ. Physiol. 282: H2296-H2308, 2002.
 *
 * The model is also described at: https://models.physiomeproject.org/e/5/bernus_wilders_zemlin_verschelde_panfilov_2002.cellml/view
 *
 * The bernus class defined here implements the Iionmodel interface and collects all routines
 * that involve the (time-dependent) gating variables. Static functions, that is functions that
 * do not dependent on the gating variables, are collected in bernus_functions.
 * For a single cell the Bernus model (and many other membrane models) is a differential equation
 *
 * \\( v' = -I_{\rm ion}(v, w) \\\
 *  w' = f(v, w) \\)
 *
 * where \\( v \\) is the membrane potential, \\( w \\) are the state variables of the membrane model
 * (''gating variables'') and both \\( I_{\rm ion} \\) and \\( f \\) are provided by the model.
 * In the Bernus model, \\( I_{\rm ion} \\) is the sum of nine different ion currents:
 * - Sodium current \\( i_{\rm Na} \\)
 * - Calcium current \\( i_{\rm Ca} \\)
 * - Transient outward current \\( i_{\rm to} \\)
 * - Delayed rectifier potassium current \\( i_{\rm K} \\)
 * - Inward rectifier potassium current \\( i_{\textrm{K}, 1} \\)
 * - Calcium background current \\( i_{\rm b, Ca} \\)
 * - Sodium background current \\( i_{\rm b, Na} \\)
 * - Sodium potassium pump \\( i_{\rm Na, K} \\)
 * - Sodium calcium pump \\( i_{\rm Na, Ca} \\)
 *
 * There are five state variables in the Bernus model that are modeled by ODEs:
 * - Sodium current gate \\( m \\)
 * - Sodium current gate \\( v \\)
 * - Calcium current gate \\( f \\)
 * - Transient outward current gate \\( to \\)
 * - Delayed rectifier potassium current gate \\( x \\)
 *
 * <b>A note on physical units</b>: The ionic current return by #ionforcing is in millivolt per millisecond \\( \\textrm{mV}~\\textrm{ms}^{-1} \\): The constants #g_na, #g_ca etc are in
 *
 * \\( \\textrm{nS}~\\textrm{pF}^{-1} = 10^3 \\textrm{S}~\\textrm{F}^{-1} = 10^3~\\textrm{s}^{-1} = \\textrm{ms}^{-1} \\)
 *
 * with S = Siemens and F = Farad, cf. Table 1 in Bernus et al. The membrane potential and equilibrium potential both are in millivolt.
 */
class bernus: public Iionmodel {
  
public:
  
  //! Constructor
  bernus();
  
  //! Destructor
  ~bernus();
  
  //! Initialize gating variables to their steady-state values
  //! for the Bernus model resting potential \\( V=-90.272 mV \\)
  void initialize(std::vector<double>* gates);
  
  double ionforcing(double,std::vector<double>*);
  
  int get_ngates();
  
  void get_gates_dt(double,  std::vector<double>*, std::vector<double>*);
  
  void rush_larsen_step(double, double,std::vector<double>*);
  
  //! Static factory function that instantiates a #bernus object and returns a pointer. Called by the #IionmodelFactory class.
  //! @param[out] Iionmodel* A pointer to an object of type #bernus.
  static Iionmodel * factory() {
    return new bernus();
  }
  
  //! Index of gating variable \\( m \\) in #gates
  static const int m_gate  = 0;
  
  //! Index of gating variable \\( v \\) in #gates
  static const int v_gate  = 1;
  
  //! Index of gating variable \\( f \\) in #gates
  static const int f_gate  = 2;
  
  //! Index of gating variable \\( to \\) in #gates
  static const int to_gate = 3;
  
  //! Index of gating variable \\( x \\) in #gates
  static const int x_gate  = 4;
  
  //! Object providing all the necessary functions to compute parameters that do not depend on the gating variables.
  static const bernus_functions bnf;
  
  //! @param[in] V Membrane potential in mV
  //! @param[in] gates Vector with values of gating variables
  //! @param[out] i_Na Sodium current
  double i_na(double V, std::vector<double>* gates);
  
  //! @param[in] V Membrane potential in mV
  //! @param[in] gates Vector with values of gating variables
  //! @param[out] i_Ca Calcium current
  double i_ca(double, std::vector<double>*);
  
  //! @param[in] V Membrane potential in mV
  //! @param[in] gates Vector with values of gating variables
  //! @param[out] i_to Transient outward current
  double i_to(double, std::vector<double>*);
  
  //! @param[in] V Membrane potential in mV
  //! @param[in] gates Vector with values of gating variables
  //! @param[out] i_k Delayed rectifier potassium current
  double i_k(double, std::vector<double>*);
  
  //! @param[in] V Membrane potential in mV
  //! @param[out] i_k1 Inward rectifier potassium current
  double i_k1(double);
  
  //! @param[in] V Membrane potential in mV
  //! @param[out] i_b_ca Calcium background current
  double i_b_ca(double);
  
  //! @param[in] V Membrane potential in mV
  //! @param[out] i_b_na Sodium background current
  double i_b_na(double);
  
  //! @param[in] V Membrane potential in mV
  //! @param[out] i_na_k Sodium potassium pump
  double i_na_k(double);
  
  //! @param[in] V Membrane potential in mV
  //! @param[out] i_na_ca Sodium calcium pump
  double i_na_ca(double);
  
  /*
   * Member variables
   */
  
  //! Constant \\( g_{\rm Na} \\) from Table 1 in Bernus et al.
  double static constexpr g_na   = 16.0;

  //! Constant \\( g_{\rm Ca} \\) from Table 1 in Bernus et al.
  double static constexpr g_ca   = 0.064;
  
  //! Constant \\( g_{\rm to} \\) from Table 1 in Bernus et al.
  double static constexpr g_to   = 0.4;
  
  //! Constant \\( g_{\rm K} \\) from Table 1 in Bernus et al.
  double static constexpr g_k    = 0.019;
  
  //! Constant \\( g_{\textrm{K},1} \\) from Table 1 in Bernus et al.
  double static constexpr g_k1   = 3.9;
  
  //! Constant \\( g_{\rm Na,b} \\) from Table 1 in Bernus et al.
  double static constexpr g_na_b = 0.001;
  
  //! Constant \\( g_{\rm Ca,b} \\) from Table 1 in Bernus et al.
  double static constexpr g_ca_b = 0.00085;
  
  //! Constant \\( g_{\rm Na,K} \\) from Table 1 in Bernus et al.
  double static constexpr g_nak  = 1.3;
  
  //! Constant \\( g_{\rm Na,Ca} \\) from Table 1 in Bernus et al.
  double static constexpr g_naca = 1000.0;
  
private:
  
  //! Number of gating variables in the Bernus model described by an ODE.
  size_t static const ngates = 5;
  
};

/*
 * To faciliate inlining, the ion current functions are implemented here, in the header file.
 * The keyword 'inline' allow e.g. the use of the -Winline flag for gcc to issue a warning if
 * the compiler was unable to actually inline the respective function.
 */

inline double bernus::ionforcing(double V, std::vector<double>* gates) {
  return i_na(V,gates)+i_ca(V,gates)+i_to(V,gates)+i_k(V,gates)+i_k1(V)+i_b_ca(V)+i_b_na(V)+i_na_k(V)+i_na_ca(V);
}

inline int bernus::get_ngates() {
  return bernus::ngates;
}

inline void bernus::get_gates_dt(double V, std::vector<double>* gates, std::vector<double>* gates_dt) {
  
  // See e.g. https://models.physiomeproject.org/e/5/bernus_wilders_zemlin_verschelde_panfilov_2002.cellml/view
  // for the ODEs for the gating variables; see also Bernus et al.
  (*gates_dt)[m_gate]  = bnf.alpha_m(V)*( 1.0 - (*gates)[m_gate])  - bnf.beta_m(V)*(*gates)[m_gate];
  (*gates_dt)[f_gate]  = bnf.alpha_f(V)*( 1.0 - (*gates)[f_gate])  - bnf.beta_f(V)*(*gates)[f_gate];
  (*gates_dt)[to_gate] = bnf.alpha_to(V)*(1.0 - (*gates)[to_gate]) - bnf.beta_to(V)*(*gates)[to_gate];

  (*gates_dt)[v_gate]  = (bnf.v_inf(V) - (*gates)[v_gate])/bnf.tau_v(V);
  (*gates_dt)[x_gate]  = (bnf.x_inf(V) - (*gates)[x_gate])/bnf.tau_x(V);
}

inline void bernus::rush_larsen_step(double V, double dt, std::vector<double>* gates) {
  
  double y_inf;
  double tau_y;
  
  // m-gate
  y_inf = bnf.alpha_m(V)/( bnf.alpha_m(V) + bnf.beta_m(V) );
  tau_y = 1.0/( bnf.alpha_m(V) + bnf.beta_m(V) );
  (*gates)[m_gate] *= exp(-dt/tau_y);
  (*gates)[m_gate] += (1.0 - exp(-dt/tau_y))*y_inf;
  
  // f-gate
  y_inf = bnf.alpha_f(V)/( bnf.alpha_f(V) + bnf.beta_f(V) );
  tau_y = 1.0/( bnf.alpha_f(V) + bnf.beta_f(V) );
  (*gates)[f_gate] *= exp(-dt/tau_y);
  (*gates)[f_gate] += (1.0 - exp(-dt/tau_y))*y_inf;
  
  // to-gate
  y_inf = bnf.alpha_to(V)/( bnf.alpha_to(V) + bnf.beta_to(V) );
  tau_y = 1.0/( bnf.alpha_to(V) + bnf.beta_to(V) );
  (*gates)[to_gate] *= exp(-dt/tau_y);
  (*gates)[to_gate] += (1.0 - exp(-dt/tau_y))*y_inf;
  
  // v-gate
  y_inf = bnf.v_inf(V);
  tau_y = bnf.tau_v(V);
  (*gates)[v_gate] *= exp(-dt/tau_y);
  (*gates)[v_gate] += (1.0 - exp(-dt/tau_y))*y_inf;
  
  // x-gate
  y_inf = bnf.x_inf(V);
  tau_y = bnf.tau_x(V);
  (*gates)[x_gate] *= exp(-dt/tau_y);
  (*gates)[x_gate] += (1.0 - exp(-dt/tau_y))*y_inf;
}


// Sodium current i_Na
inline double bernus::i_na(double V,std::vector<double>* gates){
  return g_na*pow((*gates)[m_gate], 3.0)*pow((*gates)[v_gate], 2.0)*(V - bnf.e_na);}

// Calcium current i_Ca
inline double bernus::i_ca(double V,std::vector<double>* gates){
  return g_ca*(bnf.d_inf(V))*(*gates)[f_gate]*(bnf.f_ca(V))*(V-bnf.e_ca);}

// Transient outward current i_to
inline double bernus::i_to(double V,std::vector<double>* gates){
  return g_to*(bnf.r_inf(V))*(*gates)[to_gate]*(V-bnf.e_to);}

// Delated rectifier potassium current i_K
inline double bernus::i_k(double V, std::vector<double>* gates){
  return g_k*pow( (*gates)[x_gate], 2.0)*(V-bnf.e_k);}

// Inward rectifier potassium current i_K1
inline double bernus::i_k1(double V){
  return g_k1*(bnf.k1_inf(V))*(V-bnf.e_k);}

// Calcium background current
inline double bernus::i_b_ca(double V){
  return g_ca_b*(V-bnf.e_ca);}

// Sodium background current
inline double bernus::i_b_na(double V){
  return g_na_b*(V - bnf.e_na);}

// Sodium potassium pump
inline double bernus::i_na_k(double V){
  return g_nak*(bnf.f_nak(V))*(bnf.f_nak_a(V));}

// Sodium calcium pump
inline double bernus::i_na_ca(double V){
  return g_naca*(bnf.f_naca(V));}

#endif // BERNUS_HPP