#include "mdre_parameters.hpp"
#include "mdre_node_equations.hpp"
#include "mdre_network.hpp"
#include "mdre_solver.hpp"
#include "mdre_time_series.hpp"


extern NodeParameters np;
extern Parameters p;


// constructor for TimeSeries class
TimeSeries::TimeSeries ()
{
  
  t = 0.0;
  
  timeValues.zeros(p.evalSteps);
  intensityValues0.zeros(p.evalSteps);
  intensityValues1.zeros(p.evalSteps);
  intensityValues2.zeros(p.evalSteps);
  intensityValues3.zeros(p.evalSteps);
  LyapunovValues.zeros(p.evalSteps);
  
  setIntegrationMethods();
  checkTimeSeriesFeasibility();
  setRotatingFrame();
  
  rng.seed(1);
  uniform = uniform_real_distribution<double>(0.0, 1.0);
  
}



// constructor for TimeSeries class with variable vector length
TimeSeries::TimeSeries (int vectorLength)
{
  
  t = 0.0;
  
  timeValues.zeros(vectorLength);
  intensityValues0.zeros(vectorLength);
  intensityValues1.zeros(vectorLength);
  intensityValues2.zeros(vectorLength);
  intensityValues3.zeros(vectorLength);
  LyapunovValues.zeros(p.evalSteps);
  
  setIntegrationMethods();
  checkTimeSeriesFeasibility();
  setRotatingFrame();
  
  rng.seed(1);
  uniform = uniform_real_distribution<double>(0.0, 1.0);
  
}



// adjusts numeric parameters to new dt, Ttrans, Teval and tau
void TimeSeries::updateNumericParameters ()
{
  
  p.sqrt_dt    = sqrt(p.dt / 2.0);
  p.nSteps     = int((p.Ttrans + p.Teval) / p.dt);
  p.transSteps = int(p.Ttrans / p.dt);
  p.evalSteps  = int(p.Teval / p.dt);
  p.n_tau      = int(p.tau / p.dt);
  
}



void TimeSeries::setIntegrationMethods ()
{
  
  deterministic_methods.clear();
  stochastic_methods.clear();
  
  deterministic_methods.push_back("explicitEuler");
  deterministic_methods.push_back("implicitEuler");
  deterministic_methods.push_back("explicitMidpoint");
  deterministic_methods.push_back("explicitHeun");
  deterministic_methods.push_back("explicitRungeKutta4");
  deterministic_methods.push_back("explicitEulerDelay");
  deterministic_methods.push_back("explicitRungeKutta4Delay");
  
  stochastic_methods.push_back("explicitEulerMaruyama");
  stochastic_methods.push_back("explicitEulerMaruyamaDelay");
  
  delay_methods.push_back("explicitEulerDelay");
  delay_methods.push_back("explicitRungeKutta4Delay");
  delay_methods.push_back("explicitEulerMaruyamaDelay");
  
  no_delay_methods.push_back("explicitEuler");
  no_delay_methods.push_back("implicitEuler");
  no_delay_methods.push_back("explicitMidpoint");
  no_delay_methods.push_back("explicitHeun");
  no_delay_methods.push_back("explicitRungeKutta4");
  no_delay_methods.push_back("explicitEulerMaruyama");
  
}



void TimeSeries::checkTimeSeriesFeasibility ()
{
  
  // check if correct time series mode is selected
  if (p.timeSeriesMode != "deterministic" 
   && p.timeSeriesMode != "stochastic") {
    cerr << "\nChoose time series mode from {deterministic, stochastic}.\n" << endl;
    exit(0);
  }
  
  // check that only one simulation of deterministic time series is computed
  if (p.timeSeriesMode == "deterministic" && p.nSims > 1) {
    cerr << "\nThere is no need for more than one simulation of deterministic time series.\n" << endl;
    exit(0);
  }
  
  // check if deterministic integration methods are used
  if (p.timeSeriesMode == "deterministic" 
   && find(stochastic_methods.begin(), stochastic_methods.end(), 
      p.integrationMethod) != stochastic_methods.end()) {
    cerr << "\nChoose deterministic numerical integration method from {explicitEuler, implicitEuler, explicitMidpoint, explicitHeun, explicitRungeKutta4, explicitEulerDelay, explicitRungeKutta4Delay}.\n" << endl;
    exit(0);
  }
  
  // check if stochastic integration methods are used
  else if (p.timeSeriesMode == "stochastic" 
        && find(deterministic_methods.begin(), deterministic_methods.end(), 
           p.integrationMethod) != deterministic_methods.end()) {
    cerr << "\nChoose stochastic numerical integration method from {explicitEulerMaruyama, explicitEulerMaruyamaDelay}.\n" << endl;
    exit(0);
  }
  
  // check for deterministic spontaneous emission
  if (p.timeSeriesMode == "deterministic" && p.E_sp_det == false 
      && np.beta != 0.0) {
    np.beta = 0.0;
    cout << "\nbeta is set to zero. Set command line option -E_sp_det 1 for beta > 0.\n" << endl;
  }
  
  // turn off deterministic spontaneous emission in case of stochastic
  // time series mode
  if (p.timeSeriesMode == "stochastic" && p.E_sp_det == true) {
    p.E_sp_det = false;
    cout << "\nE_sp_det is set to false because stochastic time series should be computed without deterministic spontaneous emission term.\n" << endl;
  }
  
  // set tau to zero or change network type
  if (p.networkType == "singleLaser" && p.tau != 0.0) {
    cerr << "\nNetwork type is set to singleLaser while tau is not zero. Either set tau to zero or set -networkType flag to singleLaserDelay.\n" << endl;
    exit(0);
  }
  
  // tau > 0 does not match numerical integration method
  if (p.tau != 0.0 && find(no_delay_methods.begin(), no_delay_methods.end(), 
      p.integrationMethod) != no_delay_methods.end()) {
    cerr << "\nChoose numerical integration method for ddes or sddes from {explicitEulerDelay, explicitRungeKutta4Delay, explicitEulerMaruyamaDelay}.\n" << endl;
    exit(0);
  }
  
  // negative delay error
  if (p.tau < 0.0) {
    cerr << "\nFeedback delay tau must be >= 0.\n" << endl;
    exit(0);
  }
  
}



// sets optical frequency of free running laser evaluated at fixed point
void TimeSeries::setRotatingFrame ()
{
  
  // store original parameter values and modes
  
  double dt                = p.dt;
  double Ttrans            = p.Ttrans;
  double Teval             = p.Teval;
  int    transSteps        = p.transSteps;
  int    evalSteps         = p.evalSteps;
  int    nSteps            = p.nSteps;
  int    n_tau             = p.n_tau;
  
  string networkType       = p.networkType;
  string timeSeriesMode    = p.timeSeriesMode;
  string integrationMethod = p.integrationMethod;
  string initialCondition  = p.initialCondition;
  bool   E_sp_det          = p.E_sp_det;
  
  double K_inj             = np.K_inj;
  double DnuInj            = np.DnuInj;
  
  double K_fb              = p.K_fb;
  double phi               = p.phi;
  double tau               = p.tau;
  
  
  // set parameters for computation of rotating frame
  // injection parameters are set to zero for free-running laser simulation
  
  p.dt                = 1e-3;
  p.Ttrans            = 0.0;
  p.Teval             = 10.0;
  p.nSteps            = int((p.Ttrans + p.Teval) / p.dt);
  
  p.E_sp_det          = true;
  p.networkType       = "singleLaser";
  p.timeSeriesMode    = "deterministic";
  p.integrationMethod = "explicitRungeKutta4";
  p.initialCondition  = "small";
  
  np.K_inj            = 0.0;
  np.DnuInj           = 0.0;
  
  p.K_fb              = 0.0;
  p.phi               = 0.0;
  p.tau               = 0.0;
  
  np.E0               = 0.0;
  np.omega0           = 0.0;
  
  
  Network N_rotating_frame;
  
  for (long long int j = 1; j < p.nSteps; j++) {
    timeSeriesStep(N_rotating_frame, j);
  }
  
  
  // reset time, parameters and modes to original values
  
  t = 0.0;
  
  p.dt                = dt;
  p.Ttrans            = Ttrans;
  p.Teval             = Teval;
  p.tau               = tau;
  p.transSteps        = transSteps;
  p.evalSteps         = evalSteps;
  p.nSteps            = nSteps;
  p.n_tau             = n_tau;
  
  p.E_sp_det          = E_sp_det;
  p.networkType       = networkType;
  p.timeSeriesMode    = timeSeriesMode;
  p.integrationMethod = integrationMethod;
  p.initialCondition  = initialCondition;
  
  
  // set injection parameters
  np.K_inj  = K_inj;
  np.DnuInj = DnuInj;
  
  p.K_fb    = K_fb;
  p.phi     = phi;
  p.tau     = tau;
  
  // set electric field amplitude of free-running laser in fixed point
  np.E0 = abs(N_rotating_frame.nodes[0].E);
  
  // set optical frequency of free-running laser ( = rotating frame)
  np.omega0 = N_rotating_frame.nodes[0].delta_omega();
  
}



// returns sum of electrons
double TimeSeries::electronSum ()
{
  return (N.nodes[0].w_e + 4.0 * np.N_QD * N.nodes[0].rho_ES_e 
        + 2.0 * np.N_QD * (np.f_act * N.nodes[0].rho_GS_e_act 
                         + np.f_inact * N.nodes[0].rho_GS_e_inact));
}



// returns sum of holes
double TimeSeries::holeSum ()
{
  return (N.nodes[0].w_h + 4.0 * np.N_QD * N.nodes[0].rho_ES_h 
        + 2.0 * np.N_QD * (np.f_act * N.nodes[0].rho_GS_h_act 
                         + np.f_inact * N.nodes[0].rho_GS_h_inact));
}



// returns current charge conservation state
double TimeSeries::chargeConservation ()
{
  return (electronSum() - holeSum()) / (electronSum());
}



// returns derivative of current charge conservation state
double TimeSeries::dchargeConservation ()
{
  
  dmdre = N.f(N.X, t);
  
  drho_GS_e_act   = dmdre(1).real();
  drho_GS_h_act   = dmdre(2).real();
  drho_GS_e_inact = dmdre(3).real();
  drho_GS_h_inact = dmdre(4).real();
  drho_ES_e       = dmdre(5).real();
  drho_ES_h       = dmdre(6).real();
  dw_e            = dmdre(7).real();
  dw_h            = dmdre(8).real();
  
  return (((dw_e + 4.0 * np.N_QD * drho_ES_e + 2.0 * np.N_QD 
        * (np.f_act * drho_GS_e_act + np.f_inact * drho_GS_e_inact)) 
        - (dw_h + 4.0 * np.N_QD * drho_ES_h + 2.0 * np.N_QD 
        * (np.f_act * drho_GS_h_act + np.f_inact * drho_GS_h_inact))) 
        / electronSum());
  
}



// performs time series step
void TimeSeries::timeSeriesStep (Network& _N, long long int i)
{
  
  // methods without delay
  if (p.tau == 0.0) {
    
    // ode methods
    if (p.integrationMethod == "explicitEuler") {
      S.explicitEuler(_N, t);
      if (p.Lyapunov == true) S.explicitEulerLinearized(_N, t);
    }
    
    else if (p.integrationMethod == "implicitEuler") {
      S.implicitEuler(_N, t);
    }
    
    else if (p.integrationMethod == "explicitMidpoint") {
      S.explicitMidpoint(_N, t);
    }
    
    else if (p.integrationMethod == "explicitHeun") {
      S.explicitHeun(_N, t);
    }
    
    else if (p.integrationMethod == "explicitRungeKutta4") {
      S.explicitRungeKutta4(_N, t);
      if (p.Lyapunov == true) S.explicitRungeKutta4Linearized(_N, t);
    }
    
    // sde methods
    else if (p.integrationMethod == "explicitEulerMaruyama") {
      S.explicitEulerMaruyama(_N, t);
    }
    
    else {
      cerr << "\nAppropriate numerical integration method for odes or sdes must be selected.\n" << endl;
      exit(0);
    }
    
  }
  
  // methods with delay
  else if (p.tau > 0.0) {
    
    // circular buffers for indexing
    _N.delayedStateIndex  =  i            % (p.n_tau + 1);
    _N.interpolationIndex = (i + 1)       % (p.n_tau + 1);
    _N.currentStateIndex  = (i + p.n_tau) % (p.n_tau + 1);
    
    // dde methods
    if (p.integrationMethod == "explicitEulerDelay") {
      
      S.explicitEulerDelay(_N, t);
      
      if (p.Lyapunov == true) {
        S.explicitEulerDelayLinearized(_N, t);
      }
      
    }
    
    else if (p.integrationMethod == "explicitRungeKutta4Delay") {
      
      S.explicitRungeKutta4Delay(_N, t);
      
      if (p.Lyapunov == true) {
        S.explicitRungeKutta4DelayLinearized(_N, t);
      }
      
    }
    
    // sdde methods
    else if (p.integrationMethod == "explicitEulerMaruyamaDelay") {
      S.explicitEulerMaruyamaDelay(_N, t);
    }
    
    else {
      cerr << "\nAppropriate numerical integration method for ddes or sdde must be selected.\n" << endl;
      exit(0);
    }
    
  }
  
  t = i * p.dt;
  
}



// performs time series step
void TimeSeries::writeTimeSeries (ostream& file)
{
  
  // single laser case
  if (p.networkType == "singleLaser" || p.networkType == "singleLaserDelay") {
    
    file << t                                                  << '\t' // 0
         << norm(N.nodes[0].E)                                 << '\t' // 1
         << N.nodes[0].E.real()                                << '\t' // 2
         << N.nodes[0].E.imag()                                << '\t' // 3
         << abs(N.nodes[0].E)                                  << '\t' // 4
         << arg(N.nodes[0].E)                                  << '\t' // 5
         << N.nodes[0].rho_GS_e_act                            << '\t' // 6
         << N.nodes[0].rho_GS_h_act                            << '\t' // 7
         << N.nodes[0].rho_GS_e_inact                          << '\t' // 8
         << N.nodes[0].rho_GS_h_inact                          << '\t' // 9
         << N.nodes[0].rho_ES_e                                << '\t' // 10
         << N.nodes[0].rho_ES_h                                << '\t' // 11
         << N.nodes[0].w_e                                     << '\t' // 12
         << N.nodes[0].w_h                                     << '\t' // 13
         << electronSum()                                      << '\t' // 14
         << holeSum()                                          << '\t' // 15
         << chargeConservation()                               << '\t' // 16
         << dchargeConservation()                              << '\t' // 17
         << N.nodes[0].rho_GS_e_act + N.nodes[0].rho_GS_h_act  << '\t' // 18
         << N.nodes[0].rho_ES_e + N.nodes[0].rho_ES_h          << '\t' // 19
         << N.nodes[0].S_e_rel_in() - N.nodes[0].S_e_rel_out() << '\t' // 20
         << N.nodes[0].S_h_rel_in() - N.nodes[0].S_h_rel_out() << '\t';// 21
    
    if (p.Lyapunov == true) {
      file << N.LyapunovExponent();
    }
         
    file << endl;
    
  }
  
  // four node network case
  else if (p.networkType == "fourNodeAllToAll") {
    
    file << t             << '\t' 
         << norm(N.X(0))  << '\t' << norm(N.X(9))  << '\t' 
         << norm(N.X(18)) << '\t' << norm(N.X(27)) << '\t' 
         << endl;
    
  }
  
}



// intensity of lasers in network is perturbed
void TimeSeries::networkPerturbation ()
{
  N.X(0)  += 1e-3 * uniform(rng) * exp(1i * 2.0 * M_PI * uniform(rng));
  N.X(9)  += 1e-3 * uniform(rng) * exp(1i * 2.0 * M_PI * uniform(rng));
  N.X(18) += 1e-3 * uniform(rng) * exp(1i * 2.0 * M_PI * uniform(rng));
  N.X(27) += 1e-3 * uniform(rng) * exp(1i * 2.0 * M_PI * uniform(rng));
}



// prints parameters of last simulation to console
void TimeSeries::showParameters ()
{
  
  cout << "network type:       " << p.networkType       << endl;
  cout << "time series mode:   " << p.timeSeriesMode    << endl;
  cout << "integration method: " << p.integrationMethod << endl;
  cout << "initial condition   " << p.initialCondition  << endl;
  cout << "E_sp_det:           " << p.E_sp_det  << '\n' << endl;
  
  cout << "dt                = " << p.dt                << endl;
  cout << "Ttrans            = " << p.Ttrans            << endl;
  cout << "Teval             = " << p.Teval     << '\n' << endl;
  
  cout << "J                 = " << np.J                << endl;
  cout << "sc                = " << np.sc               << endl;
  cout << "beta              = " << np.beta     << '\n' << endl;
  
  cout << "K_inj             = " << np.K_inj            << endl;
  cout << "DnuInj            = " << np.DnuInj   << '\n' << endl;
  
  cout << "K_fb              = " << p.K_fb              << endl;
  cout << "phi               = " << p.phi               << endl;
  cout << "tau               = " << p.tau       << '\n' << endl;
  
  cout << "E0                = " << np.E0               << endl;
  cout << "omega0            = " << np.omega0           << endl;
  cout << "omega_inj         = " << 2.0 * M_PI * np.DnuInj + np.omega0 << endl;
  
}

