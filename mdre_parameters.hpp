#ifndef MDRE_PARAMETERS_H
#define MDRE_PARAMETERS_H

#include <algorithm>
#include <armadillo>
#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <string>
#include <vector>


using namespace std;


class NodeParameters
{
  
  public:
  
  NodeParameters ();
  
  
  // constants
  
  const double e    = 1.6021766208e-19;		// elementary charge in C
  const double hbar = 6.582119514e-13;		// Planck constant in meVs
  const double kB   = 8.6173303e-2;			// Boltzmann constant in meVK⁻¹
  const double m0   = 9.10938356e-31;		// electron mass in kg
  const double me   = 0.043 * m0;			// effective electron mass in kg
  const double mh   = 0.45 * m0;			// effective hole mass in kg
  
  
  // control parameters
  
  double J;					// pump current rate in ns⁻¹, default value corresponds 2 * J_th for sc=1
  double sc;				// scaling factor for in-scattering rates
  double beta;				// spontaneous emission coefficient
  double K_inj;				// optical injection coupling strength
  double E0;				// magnitude of free-running laser electric field in fixed point
  double DnuInj;			// master-slave frequency detuning in GHz
  double omega0;			// angular frequency of free-running laser in ns⁻¹
  
  
  // laser parameters
  
  int    a_L;				// number of active QW-layers
  double f_act;				// fraction of resonant (active) QDs
  double f_inact;			// fraction of off-resonant (inactive) QDs
  double g_GS;				// ground state gain coefficient in ns⁻¹
  double N_QD;				// 2D InAs QD density per QW layer, normalized to 10¹¹ cm⁻²
  double R_w_loss;			// QW loss rate in ns⁻¹
  double T_e_eq;			// quasi-equilibrium electron temperature in K
  double T_h_eq;			// quasi-equilibrium hole temperature in K
  double W_GS;				// Einstein coefficient for spontaneous recombination losses of ground state in ns⁻¹
  double W_ES;				// Einstein coefficient for spontaneous recombination losses of excited state in ns⁻¹
  
  double kappa;				// cavity loss rate in ns⁻¹
  double delta_omega_ES;	// change of instantaneous frequency (excited state) in ns⁻¹
  double delta_omega_e_QW;	// change of instantaneous frequency (electrons) in ns⁻¹
  double delta_omega_h_QW;	// change of instantaneous frequency (holes) in ns⁻¹
  double eps_ES_e;			// electron excited state confinement energy in meV
  double eps_GS_e;			// electron ground state confinement energy in meV
  double eps_GS_h;			// hole ground state confinement energy in meV
  double eps_ES_h;			// hole excited state confinement energy in meV
  
  
  // calculating constant exponential factors for scattering rates 
  // to reduce computation time
  
  double factor_S_e_rel_out;
  double factor_S_h_rel_out;
  
  
  // 2D density of states
  
  const double Dens_e = me / (M_PI * hbar * hbar);
  const double Dens_h = mh / (M_PI * hbar * hbar);
  
  
  // fit parameters for charge-carrier scattering processes
  
  double A_GS_e;			// in ns⁻¹
  double B_GS_e;			// 
  double A_GS_h;			// in ns⁻¹
  double B_GS_h;			// 
  
  double A_ES_e;			// in ns⁻¹
  double B_ES_e;			// 
  double A_ES_h;			// in ns⁻¹
  double B_ES_h;			// 
  
  double C_e;				// in ns⁻¹
  double D_e;				// 
  double C_h;				// in ns⁻¹
  double D_h;				// 
  
};


class Parameters
{
  public:
  
  Parameters ();
  
  
  // parameters for numerical integration
  
  double dt;				// time step for numerical integration in ns
  double sqrt_dt;			// standard deviation for real and imaginary parts of complex gaussian noise
  double Ito_dt;			// ... of Itô-integrals
  
  double Ttrans;			// transient time in ns
  double Teval;				// evaluation time in ns
  int    outSteps;			// every outSteps'th step is written in file
  int    nSims;				// number of simulations
  
  int    nSteps;			// number of steps
  int    transSteps;		// transient steps
  int    evalSteps;			// evaluation steps
  
  
  // method mode and type parameters
  
  string networkType;
  string timeSeriesMode;
  string integrationMethod;
  string initialCondition;
  string scanType;
  bool   E_sp_det;
  bool   Lyapunov;
  
  // network coupling parameters
  
  double       K_fb;
  double       phi;
  double       tau;
  unsigned int n_tau;
  
  
  // function switches
  
  bool timeSeries;
  bool inputOutput;
  bool nonlinearScatteringRates;
  bool scatteringProcesses;
  bool injectionLineScan;
  bool injectionBifurcationLineScan;
  bool bifurcationScan;
  bool feedbackLineScan;
  bool feedbackBifurcationLineScan;
  bool rotatingFrameChange;
  bool findSNIPER;
  bool findHomoclinic;
  bool SNIPER_change;
  bool homoclinic_change;
  bool noiseCheck;
  bool TisiVarianceInjection;
  bool TisiVarianceFeedback;
  
};


#endif

