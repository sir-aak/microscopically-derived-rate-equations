#include "mdre_parameters.hpp"


// instantiation of NodeParameters and Parameters class
NodeParameters np;
Parameters p;


// NodeParameters constructor
// all parameters for MDRE system are set by instantiation
NodeParameters::NodeParameters ()
{
  
  // control parameters
  
  J                = 17.0;			// pump current rate in ns⁻¹, default value corresponds 2 * J_th for sc=1
  sc               = 1.0;			// scaling factor for in-scattering rates
  beta             = 0.0;			// spontaneous emission coefficient
  K_inj            = 0.0;			// optical injection coupling strength
  E0               = 0.0;			// magnitude of free-running laser electric field in fixed point
  DnuInj           = 0.0;			// master-slave frequency detuning in GHz
  omega0           = 0.0;			// angular frequency of free-running laser in ns⁻¹
  
  
  // laser parameters
  
  a_L              = 15;			// number of active QW-layers
  f_act            = 0.5;			// fraction of resonant (active) QDs
  f_inact          = 1.0 - f_act;	// fraction of off-resonant (inactive) QDs
  g_GS             = 230.0;			// ground state gain coefficient in ns⁻¹
  N_QD             = 1.0;			// 2D InAs QD density per QW layer, normalized to 10¹¹ cm⁻²
  R_w_loss         = 0.54;			// QW loss rate in ns⁻¹
  T_e_eq           = 300.0;			// quasi-equilibrium electron temperature in K
  T_h_eq           = 300.0;			// quasi-equilibrium hole temperature in K
  W_GS             = 0.44;			// Einstein coefficient for spontaneous recombination losses of ground state in ns⁻¹
  W_ES             = 0.55;			// Einstein coefficient for spontaneous recombination losses of excited state in ns⁻¹
  
  kappa            = 50.0;			// cavity loss rate in ns⁻¹
  delta_omega_ES   = 125.0;			// change rate of instantaneous frequency (excited state) in ns⁻¹
  delta_omega_e_QW = 11.3;			// change rate of instantaneous frequency (electrons) in ns⁻¹
  delta_omega_h_QW = 5.5;			// change rate of instantaneous frequency (holes) in ns⁻¹
  eps_ES_e         = -14.0;			// electron excited state confinement energy in meV
  eps_GS_e         = -64.0;			// electron ground state confinement energy in meV
  eps_GS_h         = -35.0;			// hole ground state confinement energy in meV
  eps_ES_h         = -15.0;			// hole excited state confinement energy in meV
  
  
  // calculating constant exponential factors for scattering rates 
  // to reduce computation time
  
  factor_S_e_rel_out = exp((eps_GS_e - eps_ES_e) / (kB * T_e_eq));
  factor_S_h_rel_out = exp((eps_GS_h - eps_ES_h) / (kB * T_h_eq));
  
  
  // fit parameters for charge-carrier scattering processes
  
  A_GS_e = 18.5;		// in ns⁻¹
  B_GS_e = 1.9;			// 
  A_GS_h = 10.5;		// in ns⁻¹
  B_GS_h = 5.3;			// 
  
  A_ES_e = 48.3;		// in ns⁻¹
  B_ES_e = 0.48;		// 
  A_ES_h = 21.4;		// in ns⁻¹
  B_ES_h = 1.8;			// 
  
  C_e    = 1014.0;		// in ns⁻¹
  D_e    = 1.4;			// 
  C_h    = 2272.0;		// in ns⁻¹
  D_h    = 2.3;			// 
  
}


// Parameters constructor
// numerical parameters, method mode and type parameters, network
// coupling parameters and function switches are set by instantiation
Parameters::Parameters ()
{
  
  // parameters for numerical integration
  
  dt         = 1e-3;
  sqrt_dt    = sqrt(dt / 2.0);
  
  Ttrans     = 0.0;
  Teval      = 50.0;
  outSteps   = 1;
  nSims      = 1;
  
  nSteps        = int((Ttrans + Teval) / dt);
  transSteps    = int(Ttrans / dt);
  evalSteps     = int(Teval / dt);
  
  
  // method mode and type parameters
  
  networkType       = "singleLaser";
  timeSeriesMode    = "deterministic";
  integrationMethod = "explicitRungeKutta4";
  initialCondition  = "small";
  scanType          = "scan";
  E_sp_det          = false;
  Lyapunov          = false;
  
  
  // network coupling parameters
  
  K_fb  = 0.0;
  phi   = 0.0;
  tau   = 0.0;
  n_tau = int(tau / dt);
  
  
  // function switches
  
  timeSeries                   = false;
  inputOutput                  = false;
  nonlinearScatteringRates     = false;
  scatteringProcesses          = false;
  injectionLineScan            = false;
  injectionBifurcationLineScan = false;
  bifurcationScan              = false;
  feedbackLineScan             = false;
  feedbackBifurcationLineScan  = false;
  rotatingFrameChange          = false;
  findSNIPER                   = false;
  findHomoclinic               = false;
  SNIPER_change                = false;
  homoclinic_change            = false;
  noiseCheck                   = false;
  TisiVarianceInjection        = false;
  TisiVarianceFeedback         = false;
  
}

