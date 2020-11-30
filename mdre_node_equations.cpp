#include "mdre_parameters.hpp"
#include "mdre_network.hpp"
#include "mdre_node_equations.hpp"


// instantiation of NodeParameters and Parameters class
extern NodeParameters np;
extern Parameters p;


// Node class constructor
Node::Node (int _nodeNumber, arma::cx_vec& X_net)
{
  
  nodeNumber = _nodeNumber;
  X_node.set_size(9);
  X_node = get_X_node(X_net);
  X_node_transformed.set_size(10);
  
}



// returns node state X_node extracted from network state X_net
arma::cx_vec Node::get_X_node (arma::cx_vec& X_net)
{
  
  X_node(0) = X_net(9 * nodeNumber);
  X_node(1) = X_net(9 * nodeNumber + 1);
  X_node(2) = X_net(9 * nodeNumber + 2);
  X_node(3) = X_net(9 * nodeNumber + 3);
  X_node(4) = X_net(9 * nodeNumber + 4);
  X_node(5) = X_net(9 * nodeNumber + 5);
  X_node(6) = X_net(9 * nodeNumber + 6);
  X_node(7) = X_net(9 * nodeNumber + 7);
  X_node(8) = X_net(9 * nodeNumber + 8);
  
  E              = X_node(0);
  rho_GS_e_act   = X_node(1).real();
  rho_GS_h_act   = X_node(2).real();
  rho_GS_e_inact = X_node(3).real();
  rho_GS_h_inact = X_node(4).real();
  rho_ES_e       = X_node(5).real();
  rho_ES_h       = X_node(6).real();
  w_e            = X_node(7).real();
  w_h            = X_node(8).real();
  
  return (X_node);
  
}



// returns right hand side of laser rate equations
// corresponds node dynamics
arma::cx_vec Node::f_node (arma::cx_vec& X_net, double t)
{
  
  X_node = get_X_node(X_net);
  
  dX_node = {
    
    // complex electric field
    (g() - 1i * (delta_omega() - np.omega0) - np.kappa) * E 
    + np.K_inj * np.kappa * np.E0 * exp(-1i * 2.0 * M_PI * np.DnuInj * t), 
    
    // electron occupation probability of active QD in ground state
    -g() * norm(E) / (np.f_act * np.a_L * np.N_QD) - R_sp_GS_act() 
    + S_GS_e_cap_act() + S_e_rel_act(), 
    
    // hole occupation probability of active QD in ground state
    -g() * norm(E) / (np.f_act * np.a_L * np.N_QD) - R_sp_GS_act() 
    + S_GS_h_cap_act() + S_h_rel_act(), 
    
    // electron occupation probability of inactive QD in ground state
    -R_sp_GS_inact() + S_GS_e_cap_inact() + S_e_rel_inact(), 
    
    // hole occupation probability of inactive QD in ground state
    -R_sp_GS_inact() + S_GS_h_cap_inact() + S_h_rel_inact(), 
    
    // electron occupation probability of excited state
    -R_sp_ES() + S_ES_e_cap() - 0.5 * (np.f_act * S_e_rel_act() 
    + np.f_inact * S_e_rel_inact()), 
    
    // hole occupation probability of excited state
    -R_sp_ES() + S_ES_h_cap() - 0.5 * (np.f_act * S_h_rel_act() 
    + np.f_inact * S_h_rel_inact()), 
    
    // 2D electron density
    -2.0 * np.N_QD * (np.f_act * S_GS_e_cap_act() 
    + np.f_inact * S_GS_e_cap_inact()) 
    - 4.0 * np.N_QD * S_ES_e_cap() - r_loss_w() + np.J, 
    
    // 2D hole density
    -2.0 * np.N_QD * (np.f_act * S_GS_h_cap_act() 
    + np.f_inact * S_GS_h_cap_inact()) 
    - 4.0 * np.N_QD * S_ES_h_cap() - r_loss_w() + np.J
    
  };
  
  // if E_sp_det is true, deterministic spontaneous emission is added
  if (p.E_sp_det == true) {
    dX_node(0) += dE_sp();
  }
  
  return (dX_node);
  
}



// returns stochastic part of right hand side of MDRE model 
// as complex Armadillo vector
arma::cx_vec Node::f_node_stoch ()
{
  dX_node_stoch = {sqrt(D()), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  return (dX_node_stoch);
}



// amplitude gain in ns⁻¹
double Node::g ()
{
  return (np.g_GS * (rho_GS_e_act + rho_GS_h_act - 1.0));
}


// change of instantaneous frequency in ns⁻¹
double Node::delta_omega ()
{
  return (np.delta_omega_ES * (rho_ES_e + rho_ES_h) 
        + np.delta_omega_e_QW * w_e + np.delta_omega_h_QW * w_h);
}


// spontaneous emission rate in ns⁻¹
double Node::D ()
{
  return (np.beta * np.a_L * np.N_QD * np.f_act * R_sp_GS_act());
}


// deterministic spontaneous emission term in ns⁻¹
complex<double> Node::dE_sp ()
{
  return (D() * E / norm(E));
}


// *********************************************************************


// electron quasi Fermi level
// e * 1e12 is a compensation factor that ensures the usage of the
// correct physical units
double Node::E_F_e_eq ()
{
  return (np.kB * np.T_e_eq * log(exp((w_e * np.e * 1e12) 
        / (np.Dens_e * np.kB * np.T_e_eq)) - 1.0));
}


// hole quasi Fermi level
// e * 1e12 is a compensation factor that ensures the usage of the
// correct physical units
double Node::E_F_h_eq ()
{
  return (np.kB * np.T_h_eq * log(exp((w_h * np.e * 1e12) 
        / (np.Dens_h * np.kB * np.T_h_eq)) - 1.0));
}


// *********************************************************************


// spontaneous recombination losses in ns⁻¹

double Node::R_sp_GS_act ()
{
  return (np.W_GS * rho_GS_e_act * rho_GS_h_act);
}


double Node::R_sp_GS_inact ()
{
  return (np.W_GS * rho_GS_e_inact * rho_GS_h_inact);
}


double Node::R_sp_ES ()
{
  return (np.W_ES * rho_ES_e * rho_ES_h);
}


// effective QW loss rate in ns⁻¹
double Node::r_loss_w ()
{
  return (np.R_w_loss * w_e * w_h);
}


// *********************************************************************


// nonlinear scattering rates in ns⁻¹:

// in-scattering rates derived from full microscopically calculated rates

// capture in-rates
double Node::S_GS_e_cap_in ()
{
  return (np.sc * np.A_GS_e * w_e * w_e / (np.B_GS_e + w_e));
}


double Node::S_GS_h_cap_in ()
{
  return (np.sc * np.A_GS_h * w_h * w_h / (np.B_GS_h + w_h));
}


double Node::S_ES_e_cap_in ()
{
  return (np.sc * np.A_ES_e * w_e * w_e / (np.B_ES_e + w_e));
}


double Node::S_ES_h_cap_in ()
{
  return (np.sc * np.A_ES_h * w_h * w_h / (np.B_ES_h + w_h));
}


// relaxation in-rates
double Node::S_e_rel_in ()
{
  return (np.sc * np.C_e * w_e / (np.D_e + w_e));
}


double Node::S_h_rel_in ()
{
  return (np.sc * np.C_h * w_h / (np.D_h + w_h));
}


// *********************************************************************


// out-scattering rates given by detailed balance relationships

// capture out-rates
double Node::S_GS_e_cap_out ()
{
  if (w_e == 0.0) return (0.0);
  
  else {
    return (S_GS_e_cap_in() 
          * exp((np.eps_GS_e - E_F_e_eq()) / (np.kB * np.T_e_eq)));
  }
}


double Node::S_GS_h_cap_out ()
{
  if (w_h == 0.0) return (0.0);
  
  else {
    return (S_GS_h_cap_in() 
          * exp((np.eps_GS_h - E_F_h_eq()) / (np.kB * np.T_h_eq)));
  }
}


double Node::S_ES_e_cap_out ()
{
  if (w_e == 0.0) return (0.0);
  
  else {
  return (S_ES_e_cap_in() 
        * exp((np.eps_ES_e - E_F_e_eq()) / (np.kB * np.T_e_eq)));
  }
}


double Node::S_ES_h_cap_out ()
{
  if (w_h == 0.0) return (0.0);
  
  else {
    return (S_ES_h_cap_in() 
          * exp((np.eps_ES_h - E_F_h_eq()) / (np.kB * np.T_h_eq)));
  }
}


// relaxation out-rates
double Node::S_e_rel_out ()
{
  return (S_e_rel_in() * np.factor_S_e_rel_out);
}


double Node::S_h_rel_out ()
{
  return (S_h_rel_in() * np.factor_S_h_rel_out);
}


// *********************************************************************


// scattering processes between different charge-carrier states in ns⁻¹:

// capture rates
double Node::S_GS_e_cap_act ()
{
  return (S_GS_e_cap_in() * (1.0 - rho_GS_e_act) 
        - S_GS_e_cap_out() * rho_GS_e_act);
}


double Node::S_GS_h_cap_act ()
{
  return (S_GS_h_cap_in() * (1.0 - rho_GS_h_act) 
        - S_GS_h_cap_out() * rho_GS_h_act);
}


double Node::S_GS_e_cap_inact ()
{
  return (S_GS_e_cap_in() * (1.0 - rho_GS_e_inact)
        - S_GS_e_cap_out() * rho_GS_e_inact);
}


double Node::S_GS_h_cap_inact ()
{
  return (S_GS_h_cap_in() * (1.0 - rho_GS_h_inact)
        - S_GS_h_cap_out() * rho_GS_h_inact);
}


double Node::S_ES_e_cap ()
{
  return (S_ES_e_cap_in() * (1.0 - rho_ES_e)
        - S_ES_e_cap_out() * rho_ES_e);
}


double Node::S_ES_h_cap ()
{
  return (S_ES_h_cap_in() * (1.0 - rho_ES_h)
        - S_ES_h_cap_out() * rho_ES_h);
}


// relaxation rates
double Node::S_e_rel_act ()
{
  return (S_e_rel_in() * rho_ES_e * (1.0 - rho_GS_e_act) 
        - S_e_rel_out() * (1.0 - rho_ES_e) * rho_GS_e_act);
}


double Node::S_h_rel_act ()
{
  return (S_h_rel_in() * rho_ES_h * (1.0 - rho_GS_h_act) 
        - S_h_rel_out() * (1.0 - rho_ES_h) * rho_GS_h_act);
}


double Node::S_e_rel_inact ()
{
  return (S_e_rel_in() * rho_ES_e * (1.0 - rho_GS_e_inact) 
        - S_e_rel_out() * (1.0 - rho_ES_e) * rho_GS_e_inact);
}


double Node::S_h_rel_inact ()
{
  return (S_h_rel_in() * rho_ES_h * (1.0 - rho_GS_h_inact) 
        - S_h_rel_out() * (1.0 - rho_ES_h) * rho_GS_h_inact);
}



// *********************************************************************
// methods for linearization
// *********************************************************************


// returns transformed node state X_node_transformed extracted from 
// transformed network state X_net_transformed
arma::colvec Node::get_X_node_transformed (arma::colvec& X_net_transformed)
{
  
  X_node_transformed(0) = X_net_transformed(10 * nodeNumber);
  X_node_transformed(1) = X_net_transformed(10 * nodeNumber + 1);
  X_node_transformed(2) = X_net_transformed(10 * nodeNumber + 2);
  X_node_transformed(3) = X_net_transformed(10 * nodeNumber + 3);
  X_node_transformed(4) = X_net_transformed(10 * nodeNumber + 4);
  X_node_transformed(5) = X_net_transformed(10 * nodeNumber + 5);
  X_node_transformed(6) = X_net_transformed(10 * nodeNumber + 6);
  X_node_transformed(7) = X_net_transformed(10 * nodeNumber + 7);
  X_node_transformed(8) = X_net_transformed(10 * nodeNumber + 8);
  X_node_transformed(9) = X_net_transformed(10 * nodeNumber + 9);
  
  d_Ex             = X_node_transformed(0);
  d_Ey             = X_node_transformed(1);
  d_rho_GS_e_act   = X_node_transformed(2);
  d_rho_GS_h_act   = X_node_transformed(3);
  d_rho_GS_e_inact = X_node_transformed(4);
  d_rho_GS_h_inact = X_node_transformed(5);
  d_rho_ES_e       = X_node_transformed(6);
  d_rho_ES_h       = X_node_transformed(7);
  d_w_e            = X_node_transformed(8);
  d_w_h            = X_node_transformed(9);
  
  return (X_node_transformed);
  
}



// returns extended real node state vector with real electric field 
// amplitude and phase instead of complex electric field
arma::colvec Node::f_node_transformed (arma::colvec& X_net_transformed, 
                                       double t)
{
  
  X_node_transformed = get_X_node_transformed(X_net_transformed);
  
  dX_node_transformed = {
    
    // real part of electric field
    (d_g() - np.kappa) * d_Ex + (d_delta_omega() - np.omega0) * d_Ey
    + np.K_inj * np.kappa * np.E0 * cos(2.0 * M_PI * np.DnuInj * t), 
    
    // imaginary part of electric field
    (d_g() - np.kappa) * d_Ey - (d_delta_omega() - np.omega0) * d_Ex
    - np.K_inj * np.kappa * np.E0 * sin(2.0 * M_PI * np.DnuInj * t), 
    
    // electron occupation probability of active QD in ground state
    -d_g() * (d_Ex * d_Ex + d_Ey * d_Ey) / (np.f_act * np.a_L * np.N_QD) 
    - d_R_sp_GS_act() + d_S_GS_e_cap_act() + d_S_e_rel_act(), 
    
    // hole occupation probability of active QD in ground state
    -d_g() * (d_Ex * d_Ex + d_Ey * d_Ey) / (np.f_act * np.a_L * np.N_QD) 
    - d_R_sp_GS_act() + d_S_GS_h_cap_act() + d_S_h_rel_act(), 
    
    // electron occupation probability of inactive QD in ground state
    -d_R_sp_GS_inact() + d_S_GS_e_cap_inact() + d_S_e_rel_inact(), 
    
    // hole occupation probability of inactive QD in ground state
    -d_R_sp_GS_inact() + d_S_GS_h_cap_inact() + d_S_h_rel_inact(), 
    
    // electron occupation probability of excited state
    -d_R_sp_ES() + d_S_ES_e_cap() - 0.5 * (np.f_act * d_S_e_rel_act() 
    + np.f_inact * d_S_e_rel_inact()), 
    
    // hole occupation probability of excited state
    -d_R_sp_ES() + d_S_ES_h_cap() - 0.5 * (np.f_act * d_S_h_rel_act() 
    + np.f_inact * d_S_h_rel_inact()), 
    
    // 2D electron density
    -2.0 * np.N_QD * (np.f_act * d_S_GS_e_cap_act() 
    + np.f_inact * d_S_GS_e_cap_inact()) 
    - 4.0 * np.N_QD * d_S_ES_e_cap() - d_r_loss_w() + np.J, 
    
    // 2D hole density
    -2.0 * np.N_QD * (np.f_act * d_S_GS_h_cap_act() 
    + np.f_inact * d_S_GS_h_cap_inact()) 
    - 4.0 * np.N_QD * d_S_ES_h_cap() - d_r_loss_w() + np.J
    
  };
  
  // if E_sp_det is true, deterministic spontaneous emission is added
  if (p.E_sp_det == true) {
    dX_node_transformed(0) += d_D() * d_Ex / (d_Ex * d_Ex + d_Ey * d_Ey);
    dX_node_transformed(1) += d_D() * d_Ey / (d_Ex * d_Ex + d_Ey * d_Ey);
  }
  
  return (dX_node_transformed);
  
}



// amplitude gain in ns⁻¹
double Node::d_g ()
{
  return (np.g_GS * (d_rho_GS_e_act + d_rho_GS_h_act - 1.0));
}


// change of instantaneous frequency in ns⁻¹
double Node::d_delta_omega ()
{
  return (np.delta_omega_ES * (d_rho_ES_e + d_rho_ES_h) 
        + np.delta_omega_e_QW * d_w_e + np.delta_omega_h_QW * d_w_h);
}


// spontaneous emission rate in ns⁻¹
double Node::d_D ()
{
  return (np.beta * np.a_L * np.N_QD * np.f_act * d_R_sp_GS_act());
}


// *********************************************************************


// electron quasi Fermi level
// e * 1e12 is a compensation factor that ensures the usage of the
// correct physical units
double Node::d_E_F_e_eq ()
{
  return (np.kB * np.T_e_eq * log(exp((d_w_e * np.e * 1e12) 
        / (np.Dens_e * np.kB * np.T_e_eq)) - 1.0));
}


// hole quasi Fermi level
// e * 1e12 is a compensation factor that ensures the usage of the
// correct physical units
double Node::d_E_F_h_eq ()
{
  return (np.kB * np.T_h_eq * log(exp((d_w_h * np.e * 1e12) 
        / (np.Dens_h * np.kB * np.T_h_eq)) - 1.0));
}


// *********************************************************************


// spontaneous recombination losses in ns⁻¹

double Node::d_R_sp_GS_act ()
{
  return (np.W_GS * d_rho_GS_e_act * d_rho_GS_h_act);
}


double Node::d_R_sp_GS_inact ()
{
  return (np.W_GS * d_rho_GS_e_inact * d_rho_GS_h_inact);
}


double Node::d_R_sp_ES ()
{
  return (np.W_ES * d_rho_ES_e * d_rho_ES_h);
}


// effective QW loss rate in ns⁻¹
double Node::d_r_loss_w ()
{
  return (np.R_w_loss * d_w_e * d_w_h);
}


// *********************************************************************


// nonlinear scattering rates in ns⁻¹:

// in-scattering rates derived from full microscopically calculated rates

// capture in-rates
double Node::d_S_GS_e_cap_in ()
{
  return (np.sc * np.A_GS_e * d_w_e * d_w_e / (np.B_GS_e + d_w_e));
}


double Node::d_S_GS_h_cap_in ()
{
  return (np.sc * np.A_GS_h * d_w_h * d_w_h / (np.B_GS_h + d_w_h));
}


double Node::d_S_ES_e_cap_in ()
{
  return (np.sc * np.A_ES_e * d_w_e * d_w_e / (np.B_ES_e + d_w_e));
}


double Node::d_S_ES_h_cap_in ()
{
  return (np.sc * np.A_ES_h * d_w_h * d_w_h / (np.B_ES_h + d_w_h));
}


// relaxation in-rates
double Node::d_S_e_rel_in ()
{
  return (np.sc * np.C_e * d_w_e / (np.D_e + d_w_e));
}


double Node::d_S_h_rel_in ()
{
  return (np.sc * np.C_h * d_w_h / (np.D_h + d_w_h));
}


// *********************************************************************


// out-scattering rates given by detailed balance relationships

// capture out-rates
double Node::d_S_GS_e_cap_out ()
{
  if (d_w_e == 0.0) return (0.0);
  
  else {
    return (d_S_GS_e_cap_in() 
          * exp((np.eps_GS_e - d_E_F_e_eq()) / (np.kB * np.T_e_eq)));
  }
}


double Node::d_S_GS_h_cap_out ()
{
  if (d_w_h == 0.0) return (0.0);
  
  else {
    return (d_S_GS_h_cap_in() 
          * exp((np.eps_GS_h - d_E_F_h_eq()) / (np.kB * np.T_h_eq)));
  }
}


double Node::d_S_ES_e_cap_out ()
{
  if (d_w_e == 0.0) return (0.0);
  
  else {
  return (d_S_ES_e_cap_in() 
        * exp((np.eps_ES_e - d_E_F_e_eq()) / (np.kB * np.T_e_eq)));
  }
}


double Node::d_S_ES_h_cap_out ()
{
  if (d_w_h == 0.0) return (0.0);
  
  else {
    return (d_S_ES_h_cap_in() 
          * exp((np.eps_ES_h - d_E_F_h_eq()) / (np.kB * np.T_h_eq)));
  }
}


// relaxation out-rates
double Node::d_S_e_rel_out ()
{
  return (d_S_e_rel_in() * np.factor_S_e_rel_out);
}


double Node::d_S_h_rel_out ()
{
  return (d_S_h_rel_in() * np.factor_S_h_rel_out);
}


// *********************************************************************


// scattering processes between different charge-carrier states in ns⁻¹:

// capture rates
double Node::d_S_GS_e_cap_act ()
{
  return (d_S_GS_e_cap_in() * (1.0 - d_rho_GS_e_act) 
        - d_S_GS_e_cap_out() * d_rho_GS_e_act);
}


double Node::d_S_GS_h_cap_act ()
{
  return (d_S_GS_h_cap_in() * (1.0 - d_rho_GS_h_act) 
        - d_S_GS_h_cap_out() * d_rho_GS_h_act);
}


double Node::d_S_GS_e_cap_inact ()
{
  return (d_S_GS_e_cap_in() * (1.0 - d_rho_GS_e_inact)
        - d_S_GS_e_cap_out() * d_rho_GS_e_inact);
}


double Node::d_S_GS_h_cap_inact ()
{
  return (d_S_GS_h_cap_in() * (1.0 - d_rho_GS_h_inact)
        - d_S_GS_h_cap_out() * d_rho_GS_h_inact);
}


double Node::d_S_ES_e_cap ()
{
  return (d_S_ES_e_cap_in() * (1.0 - d_rho_ES_e)
        - d_S_ES_e_cap_out() * d_rho_ES_e);
}


double Node::d_S_ES_h_cap ()
{
  return (d_S_ES_h_cap_in() * (1.0 - d_rho_ES_h)
        - d_S_ES_h_cap_out() * d_rho_ES_h);
}


// relaxation rates
double Node::d_S_e_rel_act ()
{
  return (d_S_e_rel_in() * d_rho_ES_e * (1.0 - d_rho_GS_e_act) 
        - d_S_e_rel_out() * (1.0 - d_rho_ES_e) * d_rho_GS_e_act);
}


double Node::d_S_h_rel_act ()
{
  return (d_S_h_rel_in() * d_rho_ES_h * (1.0 - d_rho_GS_h_act) 
        - d_S_h_rel_out() * (1.0 - d_rho_ES_h) * d_rho_GS_h_act);
}


double Node::d_S_e_rel_inact ()
{
  return (d_S_e_rel_in() * d_rho_ES_e * (1.0 - d_rho_GS_e_inact) 
        - d_S_e_rel_out() * (1.0 - d_rho_ES_e) * d_rho_GS_e_inact);
}


double Node::d_S_h_rel_inact ()
{
  return (d_S_h_rel_in() * d_rho_ES_h * (1.0 - d_rho_GS_h_inact) 
        - d_S_h_rel_out() * (1.0 - d_rho_ES_h) * d_rho_GS_h_inact);
}

