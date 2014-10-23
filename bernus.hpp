#ifndef BERNUS_HPP
#define BERNUS_HPP
//Define NDEBUG to activate runtime asserts
//#define NDEBUG
#include <algorithm>
#include <vector>
#include <assert.h>
#include "Iionmodel.hpp"
#include "bernus_functions.hpp"

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
 *
 * For a single cell the Bernus model (and many other membrane models) is a differential equation
 *
 * \\( v' = -I_{\rm ion}(v, w) \\\
 *  w' = f(v, w) \\)
 *
 * where \\( v \\) is the membrane potential, \\( w \\) are the state variables of the membrane model
 * (''gating variables'') and both \\( I_{\rm ion} \\) and \\( f \\) are provided by the model.
 *
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
 * - Sodium current \\( m \\)
 * - Sodium current \\( v \\)
 * - Calcium current \\( f \\)
 * - Transient outward current \\( to \\)
 * - Delayed rectifier potassium current gate \\( x \\)
 */
class bernus: public Iionmodel {
  
public:
  
  //! Initialize gating variables to their steady-state values
  //! for the Bernus model resting potential \\( V=-90.272 mV \\)
  bernus();
  
  //! Initializes gating variables to their steady-state values for
  //! given potential V
  //! @param[in] V Membrane potential in mV
  bernus(double);
  
  //! Destructor
  ~bernus();
  
  //! Evaluates all nine ion currents, sums them up and returns \\( I_{\rm ion} \\).
  //! @param[in] V Membrane potential in mV
  double ionforcing(double);
  
  //! Computes the time-derivative \\( f(v, w) \\) of the gating variables for
  //! the current values of \\( w \\) and a given membrane potential \\( v \\).
  //! New values are stored in member variable gates_dt
  //! @param[in] v Membrane potential in mV
  void update_gates_dt(double);
  
  //! The number of gating variables described by an ODE
  size_t static const ngates = 5;
  
  //! Gating variables \\( w \\).
  std::vector<double> gates;
  
  //! Time derivative \\( f(v,w) \\) of gating variables.
  std::vector<double> gates_dt;
  
  //! Fixed indices of the different gating variables
  static const int m_gate  = 0;
  static const int v_gate  = 1;
  static const int f_gate  = 2;
  static const int to_gate = 3;
  static const int x_gate  = 4;
  
  //! Object providing all the necessary functions to compute parameter
  static const bernus_functions bnf;
  
  /*
   * The total ionic current of the Bernus model is comprised of nine
   * different currents. In contrast to the parameter functions, the
   * ion current functions depend on the state of the cell, represented
   * by the ODE-based gating variables.
   */
  
  //! @param[in] V Membrane potential in mV
  //! @param[out] i_Na Sodium current
  double i_na(double);
  
  //! @param[in] V Membrane potential in mV
  //! @param[out] i_Ca Calcium current
  double i_ca(double);
  
  //! @param[in] V Membrane potential in mV
  //! @param[out] i_to Transient outward current
  double i_to(double);
  
  //! @param[in] V Membrane potential in mV
  //! @param[out] i_k Delayed rectifier potassium current
  double i_k(double);
  
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
  
  /**
   * Constants below are from Table 1 in Bernus et al.
   * TODO: Check physical units...
   */
  double static constexpr _g_na   = 16.0;
  double static constexpr _g_ca   = 0.064;
  double static constexpr _g_to   = 0.4;
  double static constexpr _g_k    = 0.019;
  double static constexpr _g_k1   = 3.9;
  double static constexpr _g_na_b = 0.001;
  double static constexpr _g_ca_b = 0.00085;
  double static constexpr _g_nak  = 1.3;
  double static constexpr _g_naca = 1000.0;
  
};

/*
 * To faciliate inlining, the ion current functions are implemented here, in the header file.
 * The keyword 'inline' allow e.g. the use of the -Winline flag for gcc to issue a warning if
 * the compiler was unable to actually inline the respective function.
 */

// Interface functions
inline double bernus::ionforcing(double V) {
  return i_na(V)+i_ca(V)+i_to(V)+i_k(V)+i_k1(V)+i_b_ca(V)+i_b_na(V)+i_na_k(V)+i_na_ca(V);
}

void bernus::update_gates_dt(double V) {

#ifndef NDEBUG
  // If NDEBUG is defined, make sure that all values in gates are between 0.0 and 1.0
  auto maxelem = std::max_element(std::begin(gates), std::end(gates));
  auto minelem = std::min_element(std::begin(gates), std::end(gates));
  assert( (0.0 <= *minelem) && (*maxelem <= 1.0) );
#endif
  
  // See e.g. https://models.physiomeproject.org/e/5/bernus_wilders_zemlin_verschelde_panfilov_2002.cellml/view
  // for the ODEs for the gating variables; see also Bernus et al.
  gates_dt[m_gate]  = bnf.alpha_m(V)*( 1.0 - gates[m_gate] ) - bnf.beta_m(V)*gates[m_gate];
  gates_dt[f_gate]  = bnf.alpha_f(V)*( 1.0 - gates[f_gate])  - bnf.beta_f(V)*gates[f_gate];
  gates_dt[to_gate] = bnf.alpha_to(V)*(1.0 - gates[to_gate]) - bnf.beta_to(V)*gates[to_gate];

  gates_dt[v_gate]  = (bnf.v_inf(V) - gates[v_gate])/bnf.tau_v(V);
  gates_dt[x_gate]  = (bnf.x_inf(V) - gates[x_gate])/bnf.tau_x(V);
}

// Sodium current i_Na
inline double bernus::i_na(double V){
  return _g_na*pow(gates[m_gate], 3.0)*pow(gates[v_gate], 2.0)*(V - bnf.e_na);}

// Calcium current i_Ca
inline double bernus::i_ca(double V){
  return _g_ca*(bnf.d_inf(V))*gates[f_gate]*(bnf.f_ca(V))*(V-bnf.e_ca);}

// Transient outward current i_to
inline double bernus::i_to(double V){
  return _g_to*(bnf.r_inf(V))*gates[to_gate]*(V-bnf.e_to);}

// Delated rectifier potassium current i_K
inline double bernus::i_k(double V){
  return _g_k*pow(gates[x_gate], 2.0)*(V-bnf.e_k);}

// Inward rectifier potassium current i_K1
inline double bernus::i_k1(double V){
  return _g_k1*(bnf.k1_inf(V))*(V-bnf.e_k);}

// Calcium background current
inline double bernus::i_b_ca(double V){
  return _g_ca_b*(V-bnf.e_ca);}

// Sodium background current
inline double bernus::i_b_na(double V){
  return _g_na_b*(V - bnf.e_na);}

// Sodium potassium pump
inline double bernus::i_na_k(double V){
  return _g_nak*(bnf.f_nak(V))*(bnf.f_nak_a(V));}

// Sodium calcium pump
inline double bernus::i_na_ca(double V){
  return _g_naca*(bnf.f_naca(V));}

#endif // BERNUS_HPP