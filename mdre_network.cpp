#include "mdre_parameters.hpp"
#include "mdre_node_equations.hpp"
#include "mdre_network.hpp"


extern NodeParameters np;
extern Parameters p;


// Network class constructor
Network::Network ()
{
  
  Eye.eye(10, 10);
  Df_analytic.set_size(10, 10);
  Df_numeric.set_size(10, 10);
  
  rng.seed(1);
  uniform = uniform_real_distribution<double>(0.0, 1.0);
  
  setNetworkProperties();
  setInitialConditions();
  if (p.Lyapunov == true) setInitialConditionsLinearized();
  checkNetworkFeasibility();
  
}



// matrix checks
void Network::checkNetworkFeasibility ()
{
  
  if (!A.is_square()) {
    cerr << "\nAdjacency matrix should be quadratic.\n" << endl;
    exit(0);
  }
  
  if (XHist.n_rows != A.n_rows * H.n_rows) {
    cerr << "\nState dimension must be equal to number of rows or columns of adjacency matrix times rows or columns of coupling scheme matrix.\n" << endl;
    exit(0);
  }
  
  if (XHist.n_cols != p.n_tau + 1) {
    cerr << "\nHistory array size needs to be initiated with appropriate size.\n" << endl;
    exit(0);
  }
  
}



// sets adjacency matrix A and coupling scheme H
void Network::setNetworkProperties ()
{
  
  arma::mat _A;
  arma::mat _H;
  arma::mat _Ht;
  
  // single laser without delay
  if (p.networkType == "singleLaser") {
    _A = {0.0};
  }
  
  // single laser with delay
  else if (p.networkType == "singleLaserDelay") {
    _A = {1.0};
  }
  
  // four node all-to-all network
  else if (p.networkType == "fourNodeAllToAll") {
    _A = {{0.0, 1.0, 1.0, 1.0}, 
          {1.0, 0.0, 1.0, 1.0}, 
          {1.0, 1.0, 0.0, 1.0}, 
          {1.0, 1.0, 1.0, 0.0}};
  }
  
  else {
    cerr << "\nChoose network type from {singleLaser, singleLaserDelay, fourNodeAllToAll}.\n" << endl;
    exit(0);
  }
  
  // coupling scheme
  _H.zeros(9, 9);
  _H(0, 0) = 1.0;
  
  // transformed coupling scheme
  _Ht.zeros(10, 10);
  _Ht(0, 0) =  cos(p.phi);
  _Ht(0, 1) = -sin(p.phi);
  _Ht(1, 0) =  sin(p.phi);
  _Ht(1, 1) =  cos(p.phi);
  
  A  = arma::sp_mat(_A);
  H  = arma::sp_mat(_H);
  AH = arma::kron(A, H);
  Ht = arma::sp_mat(_Ht);
  
  nNodes = A.n_rows;
  eig_gen(ew, ev, _A);
  
}



// sets initial condition
void Network::setInitialConditions ()
{
  
  arma::cx_vec X_init(A.n_rows * H.n_rows);
  
  // single laser
  if (nNodes == 1) {
    
    if (p.initialCondition == "small") {
      X_init = {1e-6, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-4, 1e-4};
    }
    
    else if (p.initialCondition == "random") {
      
      // random initial conditions
      for (unsigned int i = 0; i < A.n_rows * H.n_rows; i++) {
        
        if (i % H.n_rows == 0) {
          X_init(i) = uniform(rng) * exp(1i * 2.0 * M_PI * uniform(rng));
        }
        
        else {
          X_init(i) = uniform(rng);
        }
        
      }
      
    }
    
    else {
      cerr << "\nChoose initial condition from {synchron, asynchron, random}.\n" << endl;
      exit(0);
    }
    
  }
  
  // history array for four node all-to-all network
  else if (nNodes == 4) {
    
    if (p.initialCondition == "synchron") {
      X_init = {1e-6, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-4, 1e-4, 
                1e-6, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-4, 1e-4, 
                1e-6, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-4, 1e-4, 
                1e-6, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-4, 1e-4};
    }
    
    else if (p.initialCondition == "asynchron") {
      X_init = {1e-2, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-4, 1e-4, 
                2e-2, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-4, 1e-4, 
                3e-2, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-4, 1e-4, 
                4e-2, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-4, 1e-4};
    }
    
    else if (p.initialCondition == "random") {
      
      for (unsigned int i = 0; i < A.n_rows * H.n_rows; i++) {
        
        if (i % H.n_rows == 0) {
          X_init(i) = uniform(rng) * exp(1i * 2.0 * M_PI * uniform(rng));
        }
        
        else {
          X_init(i) = uniform(rng);
        }
        
      }
      
    }
    
    else {
      cerr << "\nChoose initial condition from {synchron, asynchron, random}.\n" << endl;
      exit(0);
    }
    
  }
  
  else {
    cerr << "\nAppropriate history array has to be implemented.\n" << endl;
    exit(0);
  }
  
  dX.set_size(X_init.n_elem);
  dX_stoch.set_size(X_init.n_elem);
  p0.set_size(X_init.n_elem);
  p1.set_size(X_init.n_elem);
  m0.set_size(X_init.n_elem);
  m1.set_size(X_init.n_elem);
  p05.set_size(X_init.n_elem);
  
  for (unsigned int i = 0; i < A.n_rows; i++) {
    Node n(i, X_init);
    nodes.push_back(n);
  }
  
  XHist.zeros(A.n_rows * H.n_rows, p.n_tau + 1);
  XHist.col(p.n_tau) = X_init;
  
  dXHist.zeros(A.n_rows * H.n_rows, p.n_tau + 1);
  dXHist.col(p.n_tau) = f_local(X_init, 0.0);		// only f_local if history ist zero
  
  if (p.tau == 0.0) {
    currentStateIndex = 0;
    X = X_init;
  }
  
  if (p.tau > 0.0) {
    
    currentStateIndex  = p.n_tau;
    delayedStateIndex  = 0;
    interpolationIndex = 1;
    
    X   = X_init;
    
    p0  = XHist.col(delayedStateIndex);
    p1  = XHist.col(interpolationIndex);
    m0  = dXHist.col(delayedStateIndex);
    m1  = dXHist.col(interpolationIndex);
    p05 = Hermite05();
    
  }
  
}



// displays eigenvalues, eigenvectors and diagonalization of A
void Network::showNetworkEigenSystem ()
{
  
  arma::mat _A = arma::mat(A);
  
  eig_gen(ew, ev, _A);
  
  arma::mat D = diagmat(real(ew));
  
  cout << "\neigenvalues:\n"  << real(ew) << endl;
  cout << "\neigenvectors:\n" << real(ev) << endl;
  
  cout << "\ndiagonalization A = S * D * S⁻¹:" << endl;
  cout << "S:\n"   << real(ev)      << endl;
  cout << "D:\n"   << real(D)       << endl;
  cout << "S⁻¹:\n" << real(inv(ev)) << endl;
  
  cout << "check: A = S * D * S⁻¹ = ...\n" 
       << real(ev) * real(D) * real(inv(ev)) << endl;
  
  cout << "\ndiagonalization D = S⁻¹ * A * S:" << endl;
  cout << "S⁻¹:\n" << real(inv(ev)) << endl;
  cout << "A:\n"   << real(_A)      << endl;
  cout << "S:\n"   << real(ev)      << endl;
  
  cout << "check: D = S⁻¹ * A * S = ...\n" 
       << real(inv(ev)) * real(_A) * real(ev) << endl;
  
}



// Hermite interpolation on unit interval at t = 0.5
arma::cx_vec Network::Hermite05 ()
{
  return (0.5 * (p0 + p1) + 0.125 * (m0 - m1));
}



// returns derivative of electric field amplitude
vector<double> Network::dAmplitude (double t)
{
  
  vector<double> dAList;
  
  for (int i = 0; i < nNodes; i++) {
    
    Amp  = abs(X(9 * nodes[i].nodeNumber));
    Ampd = abs(p0(9 * nodes[i].nodeNumber));
    phi  = arg(X(9 * nodes[i].nodeNumber));
    phid = arg(p0(9 * nodes[i].nodeNumber));
    
    dAmp = (nodes[i].g() - np.kappa) * Amp 
         + p.K_fb * np.kappa * Ampd * cos(phid + p.phi - phi) 
         + np.K_inj * np.kappa * np.E0 
         * cos(2.0 * M_PI * np.DnuInj * t + phi);
    
    if (p.E_sp_det == true) {
      dAmp += nodes[i].D() / Amp;
    }
    
    dAList.push_back(dAmp);
    
  }
  
  return (dAList);
  
}



// returns derivative of electric field phase
vector<double> Network::dphase (double t)
{
  
  vector<double> dphiList;
  
  for (int i = 0; i < nNodes; i++) {
    
    Amp  = abs(X(9 * nodes[i].nodeNumber));
    Ampd = abs(p0(9 * nodes[i].nodeNumber));
    phi  = arg(X(9 * nodes[i].nodeNumber));
    phid = arg(p0(9 * nodes[i].nodeNumber));
    
    dphi = -(nodes[i].delta_omega() - np.omega0)
         + p.K_fb * np.kappa * Ampd / Amp * sin(phid + p.phi - phi) 
         - np.K_inj * np.kappa * np.E0 / Amp 
         * sin(2.0 * M_PI * np.DnuInj * t + phi);
    
    dphiList.push_back(dphi);
    
  }
  
  return (dphiList);
  
}



// returns local deterministic part of right hand side of network
// corresponds local dynamics
arma::cx_vec Network::f_local (arma::cx_vec& _X, double t)
{
  
  // single laser case
  if (A.n_elem == 1) {
    dX = nodes[0].f_node(_X, t);
  }
  
  // network case
  else {
    
    for (unsigned int i = 0; i < A.n_rows; i++) {
      
      dX(9 * i)     = nodes[i].f_node(_X, t)(0);
      dX(9 * i + 1) = nodes[i].f_node(_X, t)(1);
      dX(9 * i + 2) = nodes[i].f_node(_X, t)(2);
      dX(9 * i + 3) = nodes[i].f_node(_X, t)(3);
      dX(9 * i + 4) = nodes[i].f_node(_X, t)(4);
      dX(9 * i + 5) = nodes[i].f_node(_X, t)(5);
      dX(9 * i + 6) = nodes[i].f_node(_X, t)(6);
      dX(9 * i + 7) = nodes[i].f_node(_X, t)(7);
      dX(9 * i + 8) = nodes[i].f_node(_X, t)(8);
      
    }
    
  }
  
  return (dX);
  
}



// returns global part of right hand side of network
// corresponds global contributions from all lasers
arma::cx_vec Network::f_global (arma::cx_vec& _X)
{
  
  // single laser case
  if (A.n_elem == 1 && A(0) == 0.0) {
    cerr << "\nThere are no global dynamics using a single laser without feedback.\n" << endl;
    exit(0);
  }
  
  // single laser with delay case
  else if (A.n_elem == 1 && A(0) == 1.0) {
    return (p.K_fb * np.kappa * exp(1i * p.phi) * AH * _X);
  }
  
  // network case
  else {
    return (p.K_fb * np.kappa * exp(1i * p.phi) * AH * _X / 3.0);
  }
  
} 



// returns stochastic part of right hand side of network
arma::cx_vec Network::f_stoch ()
{
  
  // single laser case
  if (A.n_elem == 1) {
    dX_stoch = nodes[0].f_node_stoch();
  }
  
  // network case
  else {
    
    for (unsigned int i = 0; i < A.n_rows; i++) {
      
      dX_stoch(9 * i)     = nodes[i].f_node_stoch()(0);
      dX_stoch(9 * i + 1) = nodes[i].f_node_stoch()(1);
      dX_stoch(9 * i + 2) = nodes[i].f_node_stoch()(2);
      dX_stoch(9 * i + 3) = nodes[i].f_node_stoch()(3);
      dX_stoch(9 * i + 4) = nodes[i].f_node_stoch()(4);
      dX_stoch(9 * i + 5) = nodes[i].f_node_stoch()(5);
      dX_stoch(9 * i + 6) = nodes[i].f_node_stoch()(6);
      dX_stoch(9 * i + 7) = nodes[i].f_node_stoch()(7);
      dX_stoch(9 * i + 8) = nodes[i].f_node_stoch()(8);
      
    }
    
  }
  
  return (dX_stoch);
  
}



// returns right hand side of network
// corresponds network dynamics without delay
arma::cx_vec Network::f (arma::cx_vec _X, double t)
{
  
  // single laser case
  if (A.n_elem == 1 && A(0) == 0.0) {
    return (f_local(_X, t));
  }
  
  // network case
  else {
    return (f_local(_X, t) + f_global(_X));
  }
  
}



// overloaded version for delay coupling
arma::cx_vec Network::f (arma::cx_vec _X, arma::cx_vec _Xd, double t)
{
  
  // single laser case
  if (A.n_elem == 1 && A(0) == 0.0) {
    return (f_local(_X, t));
  }
  
  // network case
  else {
    return (f_local(_X, t) + f_global(_Xd));
  }
  
}



// *********************************************************************
// methods for linearization
// *********************************************************************


//  sets initial condition for linearized system
void Network::setInitialConditionsLinearized ()
{
  
  deltaX = {uniform(rng), uniform(rng), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  
  deltap0.set_size(10);
  deltap1.set_size(10);
  deltam0.set_size(10);
  deltam1.set_size(10);
  deltap05.set_size(10);
  
  deltaXHist.zeros(10, p.n_tau + 1);
  ddeltaXHist.zeros(10, p.n_tau + 1);
  
  deltaXHist.col(currentStateIndex)  = deltaX;
  ddeltaXHist.col(currentStateIndex) = nodes[0].f_node_transformed(deltaX, 0.0);
  
  deltap0 = deltaXHist.col(delayedStateIndex);
  deltap1 = deltaXHist.col(interpolationIndex);
  deltam0 = ddeltaXHist.col(delayedStateIndex);
  deltam1 = ddeltaXHist.col(interpolationIndex);
  
}



// Hermite interpolation on unit interval at t = 0.5 for linearized system
arma::colvec Network::Hermite05Linearized ()
{
  return (0.5 * (deltap0 + deltap1) + 0.125 * (deltam0 - deltam1));
}



// sets real state vector with real and imaginary part of electric field
void Network::setX_transformed ()
{
  Xt = {X(0).real(), X(0).imag(), X(1).real(), X(2).real(), X(3).real(), 
        X(4).real(), X(5).real(), X(6).real(), X(7).real(), X(8).real()};
}



// returns Jacobian matrix of MDRE with self-feedback, optical injection
// and spontaneous emission
arma::mat Network::JacobianAnalytic (double t)
{
  
  Df_analytic = { 
    {-np.kappa - np.a_L * np.f_act * np.N_QD * np.W_GS * np.beta * X(1).real() * X(2).real() / norm(X(0)) + np.g_GS * (X(1).real() + X(2).real() - 1.0), -np.E0 * np.K_inj * np.kappa * sin(2.0 * M_PI * np.DnuInj * t + arg(X(0))) + abs(p0(0)) * p.K_fb * np.kappa * sin(p.phi - arg(X(0)) + arg(p0(0))), abs(X(0)) * np.g_GS + (np.a_L * np.f_act * np.N_QD * np.W_GS * np.beta * X(2).real()) / abs(X(0)), abs(X(0)) * np.g_GS + (np.a_L * np.f_act * np.N_QD * np.W_GS * np.beta * X(1).real()) / abs(X(0)), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
    {(np.kappa * (np.E0 * np.K_inj * sin(2.0 * M_PI * t * np.DnuInj + arg(X(0))) - abs(p0(0)) * p.K_fb * sin(p.phi - arg(X(0)) + arg(p0(0))))) / norm(X(0)), -((np.kappa * (np.E0 * np.K_inj * cos(2.0 * M_PI * t * np.DnuInj + arg(X(0))) + abs(p0(0)) * p.K_fb * cos(p.phi - arg(X(0)) + arg(p0(0))))) / abs(X(0))), 0.0, 0.0, 0.0, 0.0, -np.delta_omega_ES, -np.delta_omega_ES, -np.delta_omega_e_QW, -np.delta_omega_h_QW}, 
    {-((2.0 * abs(X(0)) * np.g_GS * (X(1).real() + X(2).real() - 1.0)) / (np.a_L * np.f_act * np.N_QD)), 0.0, -(norm(X(0)) * np.g_GS) / (np.a_L * np.f_act * np.N_QD) + np.A_GS_e * (- 1.0 - exp(np.eps_GS_e / (np.kB * np.T_e_eq)) / (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) * np.sc * X(7).real() * X(7).real()) / (np.B_GS_e + X(7).real()) + np.C_e * np.sc * X(7).real() * exp((-np.eps_ES_e + np.eps_GS_e) / (np.kB * np.T_e_eq)) * (-1.0 + X(5).real()) - X(5).real() / (np.D_e + X(7).real()) - np.W_GS * X(2).real(), -(norm(X(0)) * np.g_GS) / (np.a_L * np.f_act * np.N_QD) - np.W_GS * X(1).real(), 0.0, 0.0, np.C_e * np.sc * X(7).real() * (1.0 + (-1.0 + exp((np.eps_GS_e - np.eps_ES_e) / (np.kB * np.T_e_eq))) * X(1).real()) / (np.D_e + X(7).real()), 0.0, np.sc * ((np.A_GS_e * X(7).real() * (-1.0 * X(7).real() + 2.0 * (np.B_GS_e + X(7).real()) + X(7).real() * X(1).real() + exp(np.eps_GS_e / (np.kB * np.T_e_eq)) * X(7).real() * X(1).real() / (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) - 2.0 * (np.B_GS_e + X(7).real()) * X(1).real() - (2.0 * exp(np.eps_GS_e / (np.kB * np.T_e_eq)) * (np.B_GS_e + X(7).real()) * X(1).real()) / (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) + (1e12 * np.e * exp((1e12 * np.e * X(7).real() + np.Dens_e * np.eps_GS_e) / (np.Dens_e * np.kB * np.T_e_eq)) * X(7).real() * (np.B_GS_e + X(7).real()) * X(1).real()) / (np.Dens_e * (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) * (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) * np.kB * np.T_e_eq))) / ((np.B_GS_e + X(7).real()) * (np.B_GS_e + X(7).real())) + (np.C_e * np.D_e * exp(-np.eps_ES_e / (np.kB * np.T_e_eq)) * (exp(np.eps_ES_e / (np.kB * np.T_e_eq)) * X(5).real() * (X(1).real() - 1.0) + exp(np.eps_GS_e / (np.kB * np.T_e_eq)) * (X(5).real() - 1.0) * X(1).real())) / ((np.D_e + X(7).real()) * (np.D_e + X(7).real()))), 0.0}, 
    {-((2.0 * abs(X(0)) * np.g_GS * (X(1).real() + X(2).real() - 1.0)) / (np.a_L * np.f_act * np.N_QD)), 0.0, -norm(X(0)) * np.g_GS / (np.a_L * np.f_act * np.N_QD) - np.W_GS * X(2).real(), -norm(X(0)) * np.g_GS / (np.a_L * np.f_act * np.N_QD) + np.A_GS_h * (-1.0 - exp(np.eps_GS_h / (np.kB * np.T_h_eq)) / (-1.0 + exp((1e12 * np.e * X(8).real()) / (np.Dens_h * np.kB * np.T_h_eq)))) * np.sc * X(8).real() * X(8).real() / (np.B_GS_h + X(8).real()) + np.C_h * np.sc * X(8).real() * (exp((-np.eps_ES_h + np.eps_GS_h) / (np.kB * np.T_h_eq)) * (-1.0 + X(6).real()) - X(6).real()) / (np.D_h + X(8).real()) - np.W_GS * X(1).real(), 0.0, 0.0, 0.0, np.C_h * np.sc * X(8).real() * (1.0 + (-1.0 + exp((-np.eps_ES_h + np.eps_GS_h) / (np.kB * np.T_h_eq))) * X(2).real()) / (np.D_h + X(8).real()), 0.0, np.sc * ((np.A_GS_h * X(8).real() * (-1.0 * X(8).real() + 2.0 * (np.B_GS_h + X(8).real()) + X(8).real() * X(2).real() + exp(np.eps_GS_h / (np.kB * np.T_h_eq)) * X(8).real() * X(2).real() / (-1.0 + exp((1e12 * np.e * X(8).real()) / (np.Dens_h * np.kB * np.T_h_eq))) - 2.0 * (np.B_GS_h + X(8).real()) * X(2).real() - 2.0 * exp(np.eps_GS_h / (np.kB * np.T_h_eq)) * (np.B_GS_h + X(8).real()) * X(2).real() / (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) + 1e12 * np.e * exp((1e12 * np.e * X(8).real() + np.Dens_h * np.eps_GS_h) / (np.Dens_h * np.kB * np.T_h_eq)) * X(8).real() * (np.B_GS_h + X(8).real()) * X(2).real() / (np.Dens_h * (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) * (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) * np.kB * np.T_h_eq))) / ((np.B_GS_h + X(8).real()) * (np.B_GS_h + X(8).real())) + np.C_h * np.D_h * exp(-np.eps_ES_h / (np.kB * np.T_h_eq)) * (exp(np.eps_ES_h / (np.kB * np.T_h_eq)) * X(6).real() * (-1.0 * X(2).real()) + exp(np.eps_GS_h / (np.kB * np.T_h_eq)) * (-1.0 + X(6).real()) * X(2).real()) / ((np.D_h + X(8).real()) * (np.D_h + X(8).real())))}, 
    {0.0, 0.0, 0.0, 0.0, np.A_GS_e * (-1.0 - exp(np.eps_GS_e / (np.kB * np.T_e_eq)) / (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq)))) * np.sc * X(7).real() * X(7).real() / (np.B_GS_e + X(7).real()) + np.C_e * np.sc * X(7).real() * (exp((-np.eps_ES_e + np.eps_GS_e) / (np.kB * np.T_e_eq)) * (-1.0 + X(5).real()) - X(5).real()) / (np.D_e + X(7).real()) - np.W_GS * X(4).real(), -np.W_GS * X(3).real(), np.C_e * np.sc * X(7).real() * (1.0 +(-1.0 + exp((-np.eps_ES_e + np.eps_GS_e) / (np.kB * np.T_e_eq))) * X(3).real()) / (np.D_e + X(7).real()), 0.0, np.sc * ((np.A_GS_e * X(7).real() * (-1.0 * X(7).real() + 2.0 * (np.B_GS_e + X(7).real()) + X(7).real() * X(3).real() + exp(np.eps_GS_e / (np.kB * np.T_e_eq)) * X(7).real() * X(3).real() / (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e *np.kB * np.T_e_eq))) - 2.0 * (np.B_GS_e + X(7).real()) * X(3).real() - 2.0 * exp(np.eps_GS_e / (np.kB * np.T_e_eq)) * (np.B_GS_e + X(7).real()) * X(3).real() / (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) + 1e12 * np.e * exp((1e12 * np.e * X(7).real() + np.Dens_e * np.eps_GS_e) / (np.Dens_e * np.kB * np.T_e_eq)) * X(7).real() * (np.B_GS_e + X(7).real()) * X(3).real() / (np.Dens_e * (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) * (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) * np.kB * np.T_e_eq))) / ((np.B_GS_e + X(7).real()) * (np.B_GS_e + X(7).real())) + np.C_e * np.D_e * exp(-np.eps_ES_e / (np.kB * np.T_e_eq)) * (exp(np.eps_ES_e / (np.kB * np.T_e_eq)) * X(5).real() * (1.0 - X(3).real()) + exp(np.eps_GS_e / (np.kB * np.T_e_eq)) * (-1.0 + X(5).real()) * X(3).real()) / ((np.D_e + X(7).real()) * (np.D_e + X(7).real()))), 0.0}, 
    {0.0, 0.0, 0.0, 0.0, -np.W_GS * X(4).real(), np.A_GS_h * (-1.0 - exp(np.eps_GS_h / (np.kB * np.T_h_eq)) / (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq)))) * np.sc * X(8).real() * X(8).real() / (np.B_GS_h + X(8).real()) + np.C_h * np.sc * X(8).real() * (exp((-np.eps_ES_h + np.eps_GS_h) / (np.kB * np.T_h_eq)) * (-1.0 + X(6).real()) - X(6).real()) / (np.D_h + X(8).real()) - np.W_GS * X(3).real(), 0.0, np.C_h * np.sc * X(8).real() * (1.0 + (-1.0 + exp((-np.eps_ES_h + np.eps_GS_h) / (np.kB * np.T_h_eq))) * X(4).real()) / (np.D_h + X(8).real()), 0.0, np.sc * (np.A_GS_h * X(8).real() * (-1.0 * X(8).real() + 2.0 * (np.B_GS_h + X(8).real()) + X(8).real() * X(4).real() + exp(np.eps_GS_h / (np.kB * np.T_h_eq)) * X(8).real() * X(4).real() / (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) - 2.0 * (np.B_GS_h + X(8).real()) * X(4).real() - 2.0 * exp(np.eps_GS_h / (np.kB * np.T_h_eq)) * (np.B_GS_h + X(8).real()) * X(4).real() / (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) + 1e12 * np.e * exp((1e12 * np.e * X(8).real() + np.Dens_h * np.eps_GS_h) / (np.Dens_h * np.kB * np.T_h_eq)) * X(8).real() * (np.B_GS_h + X(8).real()) * X(4).real() / (np.Dens_h * (-1.0 + exp((1e12 * np.e * X(8).real()) / (np.Dens_h * np.kB * np.T_h_eq))) * (-1.0 + exp((1e12 * np.e * X(8).real()) / (np.Dens_h * np.kB * np.T_h_eq))) * np.kB * np.T_h_eq)) / ((np.B_GS_h + X(8).real()) * (np.B_GS_h + X(8).real())) + np.C_h * np.D_h * exp(-np.eps_ES_h / (np.kB * np.T_h_eq)) * (exp(np.eps_ES_h / (np.kB * np.T_h_eq)) * X(6).real() * (1.0 - X(4).real()) + exp(np.eps_GS_h / (np.kB * np.T_h_eq)) * (-1.0 + X(6).real()) * X(4).real()) / ((np.D_h + X(8).real()) * (np.D_h + X(8).real())))}, 
    {0.0, 0.0, -np.C_e * np.f_act * np.sc * X(7).real() * (exp((-np.eps_ES_e + np.eps_GS_e) / (np.kB * np.T_e_eq)) * (-1.0 + X(5).real()) - X(5).real()) / (2.0 * (np.D_e + X(7).real())), 0.0, -np.C_e * np.f_inact * np.sc * X(7).real() * (exp((-np.eps_ES_e + np.eps_GS_e) / (np.kB * np.T_e_eq)) * (-1.0 + X(5).real()) - X(5).real()) / (2.0 * (np.D_e + X(7).real())), 0.0, np.A_ES_e * (-1.0 - exp(np.eps_ES_e / (np.kB * np.T_e_eq)) / (-1.0 + exp((1e12 * np.e * X(7).real()) / (np.Dens_e * np.kB * np.T_e_eq)))) * np.sc * X(7).real() * X(7).real() / (np.B_ES_e + X(7).real()) - np.W_ES * X(6).real() - np.C_e * np.sc * X(7).real() * (np.f_act * (1.0 + (-1.0 + exp((-np.eps_ES_e + np.eps_GS_e) / (np.kB * np.T_e_eq))) * X(1).real()) + np.f_inact * (1.0 + (-1.0 + exp((-np.eps_ES_e + np.eps_GS_e) / (np.kB * np.T_e_eq))) * X(3).real())) / (2.0 * (np.D_e + X(7).real())), -np.W_ES * X(5).real(), 0.5 * np.sc * (np.A_ES_e * X(7).real() * (-2.0 * X(7).real() + 4.0 * (np.B_ES_e + X(7).real()) + 2.0 * X(7).real() * X(5).real() + (2.0 * exp(np.eps_ES_e / (np.kB * np.T_e_eq)) * X(7).real() * X(5).real()) / (-1.0 + exp((1e12 * np.e * X(7).real()) / (np.Dens_e * np.kB * np.T_e_eq))) - 4.0 * (np.B_ES_e + X(7).real()) * X(5).real() - 4.0 * exp(np.eps_ES_e / (np.kB * np.T_e_eq)) * (np.B_ES_e + X(7).real()) * X(5).real() / (-1.0 + exp((1e12 * np.e * X(7).real()) / (np.Dens_e * np.kB * np.T_e_eq))) + 2e12 * np.e * exp((1e12 * np.e * X(7).real() + np.Dens_e * np.eps_ES_e) / (np.Dens_e * np.kB * np.T_e_eq)) * X(7).real() * (np.B_ES_e + X(7).real()) * X(5).real() / (np.Dens_e * (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) * (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) * np.kB * np.T_e_eq)) / ((np.B_ES_e + X(7).real()) * (np.B_ES_e + X(7).real())) + np.C_e * np.D_e * exp(-(np.eps_ES_e / (np.kB * np.T_e_eq))) * (exp(np.eps_GS_e / (np.kB * np.T_e_eq)) * (1.0 - X(5).real()) * (np.f_act * X(1).real() + np.f_inact * X(3).real()) + exp(np.eps_ES_e / (np.kB * np.T_e_eq)) * X(5).real() * (np.f_act * (-1.0 + X(1).real()) + np.f_inact * (-1.0 + X(3).real()))) / ((np.D_e + X(7).real()) * (np.D_e + X(7).real()))), 0.0}, 
    {0.0, 0.0, 0.0, -np.C_h * np.f_act * np.sc * X(8).real() * (exp((-np.eps_ES_h + np.eps_GS_h) / (np.kB * np.T_h_eq)) * (-1.0 + X(6).real()) - X(6).real()) / (2.0 * (np.D_h + X(8).real())), 0.0, -np.C_h * np.f_inact * np.sc * X(8).real() * (exp((-np.eps_ES_h + np.eps_GS_h) / (np.kB * np.T_h_eq)) * (-1.0 + X(6).real()) - X(6).real()) / (2.0 * (np.D_h + X(8).real())), -np.W_ES * X(6).real(), np.A_ES_h * (-1.0 - (exp(np.eps_ES_h / (np.kB * np.T_h_eq))) / (-1.0 + exp((1e12 * np.e * X(8).real()) / (np.Dens_h * np.kB * np.T_h_eq)))) * np.sc * X(8).real() * X(8).real() / (np.B_ES_h + X(8).real()) - np.W_ES * X(5).real() - np.C_h * np.sc * X(8).real() * (np.f_act * (1.0 + (-1.0 + exp((-np.eps_ES_h + np.eps_GS_h) / (np.kB * np.T_h_eq))) * X(2).real()) + np.f_inact * (1.0 + (-1.0 + exp((-np.eps_ES_h + np.eps_GS_h) / (np.kB * np.T_h_eq))) * X(4).real())) / (2.0 * (np.D_h + X(8).real())), 0.0, 0.5 * np.sc * (np.A_ES_h * X(8).real() * (-2.0 * X(8).real() + 4.0 * (np.B_ES_h + X(8).real()) + 2.0 * X(8).real() * X(6).real() + (2.0 * exp(np.eps_ES_h / (np.kB * np.T_h_eq)) * X(8).real() * X(6).real()) / (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) - 4.0 * (np.B_ES_h + X(8).real()) * X(6).real() - 4.0 * exp(np.eps_ES_h / (np.kB * np.T_h_eq)) * (np.B_ES_h + X(8).real()) * X(6).real() / (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) + 2e12 * np.e * exp((1e12 * np.e * X(8).real() + np.Dens_h * np.eps_ES_h) / (np.Dens_h * np.kB * np.T_h_eq)) * X(8).real() * (np.B_ES_h + X(8).real()) * X(6).real() / (np.Dens_h * (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) * (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) * np.kB * np.T_h_eq)) / ((np.B_ES_h + X(8).real()) * (np.B_ES_h + X(8).real())) + np.C_h * np.D_h * exp(-np.eps_ES_h / (np.kB * np.T_h_eq)) * (exp(np.eps_GS_h / (np.kB * np.T_h_eq)) * (1.0 - X(6).real()) * (np.f_act * X(2).real() + np.f_inact * X(4).real()) + exp(np.eps_ES_h / (np.kB * np.T_h_eq)) * X(6).real() * (np.f_act * (-1.0 + X(2).real()) + np.f_inact * (-1.0 + X(4).real()))) / ((np.D_h + X(8).real()) * (np.D_h + X(8).real())))}, 
    {0.0, 0.0, -2.0 * np.A_GS_e * (-1.0 - exp(np.eps_GS_e / (np.kB * np.T_e_eq)) / (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq)))) * np.f_act * np.N_QD * np.sc * X(7).real() * X(7).real() / (np.B_GS_e + X(7).real()), 0.0, -2.0 * np.A_GS_e * (-1.0 - exp(np.eps_GS_e / (np.kB * np.T_e_eq)) / (-1.0 + exp((1e12 * np.e * X(7).real()) / (np.Dens_e * np.kB * np.T_e_eq)))) * np.f_inact * np.N_QD * np.sc * X(7).real() * X(7).real() / (np.B_GS_e + X(7).real()), 0.0, -4.0 * np.A_ES_e * (-1.0 - exp(np.eps_ES_e / (np.kB * np.T_e_eq)) / (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq)))) * np.N_QD * np.sc * X(7).real() * X(7).real() / (np.B_ES_e + X(7).real()), 0.0, -np.R_w_loss * X(8).real() + np.A_ES_e * np.N_QD * np.sc * X(7).real() * (4.0 * X(7).real() - 8.0 * (np.B_ES_e + X(7).real()) - 4.0 * X(7).real() * X(5).real() - 4.0 * exp(np.eps_ES_e / (np.kB * np.T_e_eq)) * X(7).real() * X(5).real() / (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) + 8.0 * (np.B_ES_e + X(7).real()) * X(5).real() + (8.0 * exp(np.eps_ES_e / (np.kB * np.T_e_eq)) * (np.B_ES_e + X(7).real()) * X(5).real()) / (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) - 4e12 * np.e * exp((1e12 * np.e * X(7).real() + np.Dens_e * np.eps_ES_e) / (np.Dens_e * np.kB * np.T_e_eq)) * X(7).real() * (np.B_ES_e + X(7).real()) * X(5).real() / (np.Dens_e * (-1.0 +exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) * (-1.0 +exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) * np.kB * np.T_e_eq)) / ((np.B_ES_e + X(7).real()) * (np.B_ES_e + X(7).real())) + np.A_GS_e * np.N_QD * np.sc * X(7).real() * (2.0 * np.f_act * X(7).real() + 2.0 * np.f_inact * X(7).real() - 4.0 * np.f_act * (np.B_GS_e + X(7).real()) - 4.0 * np.f_inact * (np.B_GS_e + X(7).real()) - 2.0 * np.f_act * X(7).real() * X(1).real() - 2.0 * exp(np.eps_GS_e / (np.kB * np.T_e_eq)) * np.f_act * X(7).real() * X(1).real() / (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) + 4.0 * np.f_act * (np.B_GS_e + X(7).real()) * X(1).real() + 4.0 * exp(np.eps_GS_e / (np.kB * np.T_e_eq)) * np.f_act * (np.B_GS_e + X(7).real()) * X(1).real() / (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) - 2e12 * np.e * exp((1e12 * np.e * X(7).real() + np.Dens_e * np.eps_GS_e) / (np.Dens_e * np.kB * np.T_e_eq)) * np.f_act * X(7).real() * (np.B_GS_e + X(7).real()) * X(1).real() / (np.Dens_e * (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) * (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) * np.kB * np.T_e_eq) - 2.0 * np.f_inact * X(7).real() * X(3).real() - 2.0 * exp(np.eps_GS_e / (np.kB * np.T_e_eq)) * np.f_inact * X(7).real() * X(3).real() / (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) + 4.0 * np.f_inact * (np.B_GS_e + X(7).real()) * X(3).real() + 4.0 * exp(np.eps_GS_e / (np.kB * np.T_e_eq)) * np.f_inact * (np.B_GS_e + X(7).real()) * X(3).real() / (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) - 2e12 * np.e * exp((1e12 * np.e * X(7).real() + np.Dens_e * np.eps_GS_e) / (np.Dens_e * np.kB * np.T_e_eq)) * np.f_inact * X(7).real() * (np.B_GS_e + X(7).real()) * X(3).real() / (np.Dens_e * (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) * (-1.0 + exp(1e12 * np.e * X(7).real() / (np.Dens_e * np.kB * np.T_e_eq))) * np.kB * np.T_e_eq)) / ((np.B_GS_e + X(7).real()) * (np.B_GS_e + X(7).real())), -np.R_w_loss * X(7).real()}, 
    {0.0, 0.0, 0.0, -2.0 * np.A_GS_h * (-1.0 - exp(np.eps_GS_h / (np.kB * np.T_h_eq)) / (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq)))) * np.f_act * np.N_QD * np.sc * X(8).real() * X(8).real() / (np.B_GS_h + X(8).real()), 0.0, -2.0 * np.A_GS_h * (-1.0 - exp(np.eps_GS_h / (np.kB * np.T_h_eq)) / (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq)))) * np.f_inact * np.N_QD * np.sc * X(8).real() * X(8).real() / (np.B_GS_h + X(8).real()), 0.0, -4.0 * np.A_ES_h * (-1.0 - exp(np.eps_ES_h / (np.kB * np.T_h_eq)) / (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq)))) * np.N_QD * np.sc * X(8).real() * X(8).real() / (np.B_ES_h + X(8).real()), -np.R_w_loss * X(8).real(), -np.R_w_loss * X(7).real() + np.N_QD * np.sc * X(8).real() * (np.A_ES_h * (4.0 * X(8).real() - 8.0 * (np.B_ES_h + X(8).real()) - 4.0 * X(8).real() * X(6).real() - 4.0 * exp(np.eps_ES_h / (np.kB * np.T_h_eq)) * X(8).real() * X(6).real() / (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) + 8.0 * (np.B_ES_h + X(8).real()) * X(6).real() + 8.0 * exp(np.eps_ES_h / (np.kB * np.T_h_eq)) * (np.B_ES_h + X(8).real()) * X(6).real() / (-1.0 + exp((1e12 * np.e * X(8).real()) / (np.Dens_h * np.kB * np.T_h_eq))) - 4e12 * np.e * exp((1e12 * np.e * X(8).real() + np.Dens_h * np.eps_ES_h) / (np.Dens_h * np.kB * np.T_h_eq)) * X(8).real() * (np.B_ES_h + X(8).real()) * X(6).real() / (np.Dens_h * (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) * (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) * np.kB * np.T_h_eq)) / ((np.B_ES_h + X(8).real()) * (np.B_ES_h + X(8).real())) + np.A_GS_h * (2.0 * np.f_act * X(8).real() + 2.0 * np.f_inact * X(8).real() - 4.0 * np.f_act * (np.B_GS_h + X(8).real()) - 4.0 * np.f_inact * (np.B_GS_h + X(8).real()) - 2.0 * np.f_act * X(8).real() * X(2).real() - 2.0 * exp(np.eps_GS_h / (np.kB * np.T_h_eq)) * np.f_act * X(8).real() * X(2).real() / (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) + 4.0 * np.f_act * (np.B_GS_h + X(8).real()) * X(2).real() + 4.0 * exp(np.eps_GS_h / (np.kB * np.T_h_eq)) * np.f_act * (np.B_GS_h + X(8).real()) * X(2).real() / (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) - 2e12 * np.e * exp((1e12 * np.e * X(8).real() + np.Dens_h * np.eps_GS_h) / (np.Dens_h * np.kB * np.T_h_eq)) * np.f_act * X(8).real() * (np.B_GS_h + X(8).real()) * X(2).real() / (np.Dens_h * (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) * (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) * np.kB * np.T_h_eq) - 2.0 * np.f_inact * X(8).real() * X(4).real() - 2.0 * exp(np.eps_GS_h / (np.kB * np.T_h_eq)) * np.f_inact * X(8).real() * X(4).real() / (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) + 4.0 * np.f_inact * (np.B_GS_h + X(8).real()) * X(4).real() + 4.0 * exp(np.eps_GS_h / (np.kB * np.T_h_eq)) * np.f_inact * (np.B_GS_h + X(8).real()) * X(4).real() / (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) - 2e12 * np.e * exp((1e12 * np.e * X(8).real() + np.Dens_h * np.eps_GS_h) / (np.Dens_h * np.kB * np.T_h_eq)) * np.f_inact * X(8).real() * (np.B_GS_h + X(8).real()) * X(4).real() / (np.Dens_h * (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) * (-1.0 + exp(1e12 * np.e * X(8).real() / (np.Dens_h * np.kB * np.T_h_eq))) * np.kB * np.T_h_eq)) / ((np.B_GS_h + X(8).real()) * (np.B_GS_h + X(8).real())))}
  };
  
  return (Df_analytic);
  
}



// returns numerical Jacobian based on finite differences
arma::mat Network::JacobianNumeric (double t, string differenceType, 
                                    double h)
{
  
  setX_transformed();
  
  // forward difference quotient
  if (differenceType == "forward") {
    
    XEx              = Xt + h * Eye.col(0);
    XEy              = Xt + h * Eye.col(1);
    Xrho_GS_e_act    = Xt + h * Eye.col(2);
    Xrho_GS_h_act    = Xt + h * Eye.col(3);
    Xrho_GS_e_inact  = Xt + h * Eye.col(4);
    Xrho_GS_h_inact  = Xt + h * Eye.col(5);
    Xrho_ES_e        = Xt + h * Eye.col(6);
    Xrho_ES_h        = Xt + h * Eye.col(7);
    Xw_e             = Xt + h * Eye.col(8);
    Xw_h             = Xt + h * Eye.col(9);
    
    fX               = nodes[0].f_node_transformed(Xt, t);
    
    fXEx             = nodes[0].f_node_transformed(XEx, t);
    fXEy             = nodes[0].f_node_transformed(XEy, t);
    fXrho_GS_e_act   = nodes[0].f_node_transformed(Xrho_GS_e_act, t);
    fXrho_GS_h_act   = nodes[0].f_node_transformed(Xrho_GS_h_act, t);
    fXrho_GS_e_inact = nodes[0].f_node_transformed(Xrho_GS_e_inact, t);
    fXrho_GS_h_inact = nodes[0].f_node_transformed(Xrho_GS_h_inact, t);
    fXrho_ES_e       = nodes[0].f_node_transformed(Xrho_ES_e, t);
    fXrho_ES_h       = nodes[0].f_node_transformed(Xrho_ES_h, t);
    fXw_e            = nodes[0].f_node_transformed(Xw_e, t);
    fXw_h            = nodes[0].f_node_transformed(Xw_h, t);
    
    Df_numeric = {
      {fXEx(0) - fX(0), fXEy(0) - fX(0), fXrho_GS_e_act(0) - fX(0), fXrho_GS_h_act(0) - fX(0), fXrho_GS_e_inact(0) - fX(0), fXrho_GS_h_inact(0) - fX(0), fXrho_ES_e(0) - fX(0), fXrho_ES_h(0) - fX(0), fXw_e(0) - fX(0), fXw_h(0) - fX(0)}, 
      {fXEx(1) - fX(1), fXEy(1) - fX(1), fXrho_GS_e_act(1) - fX(1), fXrho_GS_h_act(1) - fX(1), fXrho_GS_e_inact(1) - fX(1), fXrho_GS_h_inact(1) - fX(1), fXrho_ES_e(1) - fX(1), fXrho_ES_h(1) - fX(1), fXw_e(1) - fX(1), fXw_h(1) - fX(1)}, 
      {fXEx(2) - fX(2), fXEy(2) - fX(2), fXrho_GS_e_act(2) - fX(2), fXrho_GS_h_act(2) - fX(2), fXrho_GS_e_inact(2) - fX(2), fXrho_GS_h_inact(2) - fX(2), fXrho_ES_e(2) - fX(2), fXrho_ES_h(2) - fX(2), fXw_e(2) - fX(2), fXw_h(2) - fX(2)}, 
      {fXEx(3) - fX(3), fXEy(3) - fX(3), fXrho_GS_e_act(3) - fX(3), fXrho_GS_h_act(3) - fX(3), fXrho_GS_e_inact(3) - fX(3), fXrho_GS_h_inact(3) - fX(3), fXrho_ES_e(3) - fX(3), fXrho_ES_h(3) - fX(3), fXw_e(3) - fX(3), fXw_h(3) - fX(3)}, 
      {fXEx(4) - fX(4), fXEy(4) - fX(4), fXrho_GS_e_act(4) - fX(4), fXrho_GS_h_act(4) - fX(4), fXrho_GS_e_inact(4) - fX(4), fXrho_GS_h_inact(4) - fX(4), fXrho_ES_e(4) - fX(4), fXrho_ES_h(4) - fX(4), fXw_e(4) - fX(4), fXw_h(4) - fX(4)}, 
      {fXEx(5) - fX(5), fXEy(5) - fX(5), fXrho_GS_e_act(5) - fX(5), fXrho_GS_h_act(5) - fX(5), fXrho_GS_e_inact(5) - fX(5), fXrho_GS_h_inact(5) - fX(5), fXrho_ES_e(5) - fX(5), fXrho_ES_h(5) - fX(5), fXw_e(5) - fX(5), fXw_h(5) - fX(5)}, 
      {fXEx(6) - fX(6), fXEy(6) - fX(6), fXrho_GS_e_act(6) - fX(6), fXrho_GS_h_act(6) - fX(6), fXrho_GS_e_inact(6) - fX(6), fXrho_GS_h_inact(6) - fX(6), fXrho_ES_e(6) - fX(6), fXrho_ES_h(6) - fX(6), fXw_e(6) - fX(6), fXw_h(6) - fX(6)}, 
      {fXEx(7) - fX(7), fXEy(7) - fX(7), fXrho_GS_e_act(7) - fX(7), fXrho_GS_h_act(7) - fX(7), fXrho_GS_e_inact(7) - fX(7), fXrho_GS_h_inact(7) - fX(7), fXrho_ES_e(7) - fX(7), fXrho_ES_h(7) - fX(7), fXw_e(7) - fX(7), fXw_h(7) - fX(7)}, 
      {fXEx(8) - fX(8), fXEy(8) - fX(8), fXrho_GS_e_act(8) - fX(8), fXrho_GS_h_act(8) - fX(8), fXrho_GS_e_inact(8) - fX(8), fXrho_GS_h_inact(8) - fX(8), fXrho_ES_e(8) - fX(8), fXrho_ES_h(8) - fX(8), fXw_e(8) - fX(8), fXw_h(8) - fX(8)}, 
      {fXEx(9) - fX(9), fXEy(9) - fX(9), fXrho_GS_e_act(9) - fX(9), fXrho_GS_h_act(9) - fX(9), fXrho_GS_e_inact(9) - fX(9), fXrho_GS_h_inact(9) - fX(9), fXrho_ES_e(9) - fX(9), fXrho_ES_h(9) - fX(9), fXw_e(9) - fX(9), fXw_h(9) - fX(9)}
    };
    
  }
  
  // central difference quotient
  else if (differenceType == "central") {
    
    XEx_f              = Xt + 0.5 * h * Eye.col(0);
    XEy_f              = Xt + 0.5 * h * Eye.col(1);
    Xrho_GS_e_act_f    = Xt + 0.5 * h * Eye.col(2);
    Xrho_GS_h_act_f    = Xt + 0.5 * h * Eye.col(3);
    Xrho_GS_e_inact_f  = Xt + 0.5 * h * Eye.col(4);
    Xrho_GS_h_inact_f  = Xt + 0.5 * h * Eye.col(5);
    Xrho_ES_e_f        = Xt + 0.5 * h * Eye.col(6);
    Xrho_ES_h_f        = Xt + 0.5 * h * Eye.col(7);
    Xw_e_f             = Xt + 0.5 * h * Eye.col(8);
    Xw_h_f             = Xt + 0.5 * h * Eye.col(9);
    
    XEx_b              = Xt - 0.5 * h * Eye.col(0);
    XEy_b              = Xt - 0.5 * h * Eye.col(1);
    Xrho_GS_e_act_b    = Xt - 0.5 * h * Eye.col(2);
    Xrho_GS_h_act_b    = Xt - 0.5 * h * Eye.col(3);
    Xrho_GS_e_inact_b  = Xt - 0.5 * h * Eye.col(4);
    Xrho_GS_h_inact_b  = Xt - 0.5 * h * Eye.col(5);
    Xrho_ES_e_b        = Xt - 0.5 * h * Eye.col(6);
    Xrho_ES_h_b        = Xt - 0.5 * h * Eye.col(7);
    Xw_e_b             = Xt - 0.5 * h * Eye.col(8);
    Xw_h_b             = Xt - 0.5 * h * Eye.col(9);
    
    fXEx_f             = nodes[0].f_node_transformed(XEx_f, t);
    fXEy_f             = nodes[0].f_node_transformed(XEy_f, t);
    fXrho_GS_e_act_f   = nodes[0].f_node_transformed(Xrho_GS_e_act_f, t);
    fXrho_GS_h_act_f   = nodes[0].f_node_transformed(Xrho_GS_h_act_f, t);
    fXrho_GS_e_inact_f = nodes[0].f_node_transformed(Xrho_GS_e_inact_f, t);
    fXrho_GS_h_inact_f = nodes[0].f_node_transformed(Xrho_GS_h_inact_f, t);
    fXrho_ES_e_f       = nodes[0].f_node_transformed(Xrho_ES_e_f, t);
    fXrho_ES_h_f       = nodes[0].f_node_transformed(Xrho_ES_h_f, t);
    fXw_e_f            = nodes[0].f_node_transformed(Xw_e_f, t);
    fXw_h_f            = nodes[0].f_node_transformed(Xw_h_f, t);
    
    fXEx_b             = nodes[0].f_node_transformed(XEx_b, t);
    fXEy_b             = nodes[0].f_node_transformed(XEy_b, t);
    fXrho_GS_e_act_b   = nodes[0].f_node_transformed(Xrho_GS_e_act_b, t);
    fXrho_GS_h_act_b   = nodes[0].f_node_transformed(Xrho_GS_h_act_b, t);
    fXrho_GS_e_inact_b = nodes[0].f_node_transformed(Xrho_GS_e_inact_b, t);
    fXrho_GS_h_inact_b = nodes[0].f_node_transformed(Xrho_GS_h_inact_b, t);
    fXrho_ES_e_b       = nodes[0].f_node_transformed(Xrho_ES_e_b, t);
    fXrho_ES_h_b       = nodes[0].f_node_transformed(Xrho_ES_h_b, t);
    fXw_e_b            = nodes[0].f_node_transformed(Xw_e_b, t);
    fXw_h_b            = nodes[0].f_node_transformed(Xw_h_b, t);
    
    Df_numeric = {
      {fXEx_f(0) - fXEx_b(0), fXEy_f(0) - fXEy_b(0), fXrho_GS_e_act_f(0) - fXrho_GS_e_act_b(0), fXrho_GS_h_act_f(0) - fXrho_GS_h_act_b(0), fXrho_GS_e_inact_f(0) - fXrho_GS_e_inact_b(0), fXrho_GS_h_inact_f(0) - fXrho_GS_h_inact_b(0), fXrho_ES_e_f(0) - fXrho_ES_e_b(0), fXrho_ES_h_f(0) - fXrho_ES_h_b(0), fXw_e_f(0) - fXw_e_b(0), fXw_h_f(0) - fXw_h_b(0)}, 
      {fXEx_f(1) - fXEx_b(1), fXEy_f(1) - fXEy_b(1), fXrho_GS_e_act_f(1) - fXrho_GS_e_act_b(1), fXrho_GS_h_act_f(1) - fXrho_GS_h_act_b(1), fXrho_GS_e_inact_f(1) - fXrho_GS_e_inact_b(1), fXrho_GS_h_inact_f(1) - fXrho_GS_h_inact_b(1), fXrho_ES_e_f(1) - fXrho_ES_e_b(1), fXrho_ES_h_f(1) - fXrho_ES_h_b(1), fXw_e_f(1) - fXw_e_b(1), fXw_h_f(1) - fXw_h_b(1)}, 
      {fXEx_f(2) - fXEx_b(2), fXEy_f(2) - fXEy_b(2), fXrho_GS_e_act_f(2) - fXrho_GS_e_act_b(2), fXrho_GS_h_act_f(2) - fXrho_GS_h_act_b(2), fXrho_GS_e_inact_f(2) - fXrho_GS_e_inact_b(2), fXrho_GS_h_inact_f(2) - fXrho_GS_h_inact_b(2), fXrho_ES_e_f(2) - fXrho_ES_e_b(2), fXrho_ES_h_f(2) - fXrho_ES_h_b(2), fXw_e_f(2) - fXw_e_b(2), fXw_h_f(2) - fXw_h_b(2)}, 
      {fXEx_f(3) - fXEx_b(3), fXEy_f(3) - fXEy_b(3), fXrho_GS_e_act_f(3) - fXrho_GS_e_act_b(3), fXrho_GS_h_act_f(3) - fXrho_GS_h_act_b(3), fXrho_GS_e_inact_f(3) - fXrho_GS_e_inact_b(3), fXrho_GS_h_inact_f(3) - fXrho_GS_h_inact_b(3), fXrho_ES_e_f(3) - fXrho_ES_e_b(3), fXrho_ES_h_f(3) - fXrho_ES_h_b(3), fXw_e_f(3) - fXw_e_b(3), fXw_h_f(3) - fXw_h_b(3)}, 
      {fXEx_f(4) - fXEx_b(4), fXEy_f(4) - fXEy_b(4), fXrho_GS_e_act_f(4) - fXrho_GS_e_act_b(4), fXrho_GS_h_act_f(4) - fXrho_GS_h_act_b(4), fXrho_GS_e_inact_f(4) - fXrho_GS_e_inact_b(4), fXrho_GS_h_inact_f(4) - fXrho_GS_h_inact_b(4), fXrho_ES_e_f(4) - fXrho_ES_e_b(4), fXrho_ES_h_f(4) - fXrho_ES_h_b(4), fXw_e_f(4) - fXw_e_b(4), fXw_h_f(4) - fXw_h_b(4)}, 
      {fXEx_f(5) - fXEx_b(5), fXEy_f(5) - fXEy_b(5), fXrho_GS_e_act_f(5) - fXrho_GS_e_act_b(5), fXrho_GS_h_act_f(5) - fXrho_GS_h_act_b(5), fXrho_GS_e_inact_f(5) - fXrho_GS_e_inact_b(5), fXrho_GS_h_inact_f(5) - fXrho_GS_h_inact_b(5), fXrho_ES_e_f(5) - fXrho_ES_e_b(5), fXrho_ES_h_f(5) - fXrho_ES_h_b(5), fXw_e_f(5) - fXw_e_b(5), fXw_h_f(5) - fXw_h_b(5)}, 
      {fXEx_f(6) - fXEx_b(6), fXEy_f(6) - fXEy_b(6), fXrho_GS_e_act_f(6) - fXrho_GS_e_act_b(6), fXrho_GS_h_act_f(6) - fXrho_GS_h_act_b(6), fXrho_GS_e_inact_f(6) - fXrho_GS_e_inact_b(6), fXrho_GS_h_inact_f(6) - fXrho_GS_h_inact_b(6), fXrho_ES_e_f(6) - fXrho_ES_e_b(6), fXrho_ES_h_f(6) - fXrho_ES_h_b(6), fXw_e_f(6) - fXw_e_b(6), fXw_h_f(6) - fXw_h_b(6)}, 
      {fXEx_f(7) - fXEx_b(7), fXEy_f(7) - fXEy_b(7), fXrho_GS_e_act_f(7) - fXrho_GS_e_act_b(7), fXrho_GS_h_act_f(7) - fXrho_GS_h_act_b(7), fXrho_GS_e_inact_f(7) - fXrho_GS_e_inact_b(7), fXrho_GS_h_inact_f(7) - fXrho_GS_h_inact_b(7), fXrho_ES_e_f(7) - fXrho_ES_e_b(7), fXrho_ES_h_f(7) - fXrho_ES_h_b(7), fXw_e_f(7) - fXw_e_b(7), fXw_h_f(7) - fXw_h_b(7)}, 
      {fXEx_f(8) - fXEx_b(8), fXEy_f(8) - fXEy_b(8), fXrho_GS_e_act_f(8) - fXrho_GS_e_act_b(8), fXrho_GS_h_act_f(8) - fXrho_GS_h_act_b(8), fXrho_GS_e_inact_f(8) - fXrho_GS_e_inact_b(8), fXrho_GS_h_inact_f(8) - fXrho_GS_h_inact_b(8), fXrho_ES_e_f(8) - fXrho_ES_e_b(8), fXrho_ES_h_f(8) - fXrho_ES_h_b(8), fXw_e_f(8) - fXw_e_b(8), fXw_h_f(8) - fXw_h_b(8)}, 
      {fXEx_f(9) - fXEx_b(9), fXEy_f(9) - fXEy_b(9), fXrho_GS_e_act_f(9) - fXrho_GS_e_act_b(9), fXrho_GS_h_act_f(9) - fXrho_GS_h_act_b(9), fXrho_GS_e_inact_f(9) - fXrho_GS_e_inact_b(9), fXrho_GS_h_inact_f(9) - fXrho_GS_h_inact_b(9), fXrho_ES_e_f(9) - fXrho_ES_e_b(9), fXrho_ES_h_f(9) - fXrho_ES_h_b(9), fXw_e_f(9) - fXw_e_b(9), fXw_h_f(9) - fXw_h_b(9)}
    };
    
  }
  
  else {
    cerr << "\nChoose finite difference type from {forward, central}.\n" << endl;
    exit(0);
  }
  
  return (Df_numeric / h);
  
}



// returns local deterministic part of right hand side of linearized
// system corresponds local dynamics
arma::colvec Network::f_local_linearized (arma::colvec& _deltaX, double t)
{
  return (JacobianNumeric(t, "central") * _deltaX);
}



// returns global part of right hand side of linearized system
// corresponds global contributions from all lasers
arma::colvec Network::f_global_linearized (arma::colvec& _deltaX)
{
  return (-1.0 * p.K_fb * np.kappa * Ht * _deltaX / 3.0);
}



// returns right hand side of linearized system
// corresponds network dynamics without delay
arma::colvec Network::f_linearized (arma::colvec _deltaX, double t)
{
  return (f_local_linearized(_deltaX, t) + f_global_linearized(_deltaX));
}



// overloaded version for delay coupling
arma::colvec Network::f_linearized (arma::colvec _deltaX, 
                                    arma::colvec _deltaXd, double t)
{
  return (f_local_linearized(_deltaX, t) + f_global_linearized(_deltaXd));
}



// returns Lyapunov exponent
double Network::LyapunovExponent ()
{
  return (log(norm(deltaX)) / p.dt);
}

