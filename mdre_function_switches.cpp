#include "mdre_parameters.hpp"
#include "mdre_network.hpp"
#include "mdre_solver.hpp"
#include "mdre_auxiliary.hpp"
#include "mdre_analyzer.hpp"
#include "mdre_TISI_finder.hpp"
#include "mdre_features.hpp"
#include "mdre_function_switches.hpp"


// enables function calls / usage of MDRE features in command line
void enableFunctionSwitches (int argc, char* argv[])
{
  
  // links to instantiations of NodeParameters and Parameters class 
  // in mdre_processes.cpp
  extern NodeParameters np;
  extern Parameters p;
  
  
  // function switch for ..
  
  // time series
  if (p.timeSeries == true) {
    
    p.networkType       = "fourNodeAllToAll";
    p.timeSeriesMode    = "stochastic";
    p.integrationMethod = "explicitEulerMaruyamaDelay";
    p.initialCondition  = "synchron";
    //~ p.E_sp_det          = true;
    //~ p.Lyapunov          = true;
    
    p.E_sp_det            = false;
    p.dt                  = 1e-4;
    p.Ttrans              = 0.0;
    p.Teval               = 1000.0;
    
    p.K_fb                = 0.05;
    p.tau                 = 8e-2;
    p.phi                 = M_PI / 2.0;
    np.K_inj              = 0.05;
    
    timeSeries();
    
  }
  
  
  // input-output
  if (p.inputOutput == true) {
    
    p.networkType       = "singleLaser";
    p.timeSeriesMode    = "deterministic";
    p.integrationMethod = "explicitRungeKutta4";
    p.initialCondition  = "small";
    p.E_sp_det          = false;
    
    p.Ttrans = 0.0;
    p.Teval  = 200.0;
    
    inputOutput_mdre();
    
  }
  
  
  // nonlinear scattering rates
  if (p.nonlinearScatteringRates == true) {
    
    p.networkType       = "singleLaser";
    p.timeSeriesMode    = "deterministic";
    p.integrationMethod = "explicitRungeKutta4";
    p.initialCondition  = "small";
    p.E_sp_det          = false;
    
    nonlinearScatteringRates();
    
  }
  
  
  // scattering processes
  if (p.scatteringProcesses == true) {
    
    p.networkType       = "singleLaser";
    p.timeSeriesMode    = "deterministic";
    p.integrationMethod = "explicitRungeKutta4";
    p.initialCondition  = "small";
    p.E_sp_det          = false;
    
    scatteringProcesses();
    
  }
  
  
  // injection line-scan
  if (p.injectionLineScan == true) {
    
    p.networkType       = "singleLaserDelay";
    p.timeSeriesMode    = "deterministic";
    p.integrationMethod = "explicitRungeKutta4Delay";
    p.initialCondition  = "small";
    p.E_sp_det          = false;
    p.Lyapunov          = true;
    
    p.Ttrans = 1000.0;
    p.Teval  = 50.0;
    
    p.K_fb   = 0.05;
    p.phi    = 0.0;
    p.tau    = 8e-2;
    
    p.scanType = "scan";
    
    injectionLineScan();
    
  }
  
  
  // injection bifurcation line-scan
  if (p.injectionBifurcationLineScan == true) {
    
    p.networkType       = "singleLaserDelay";
    p.timeSeriesMode    = "deterministic";
    p.integrationMethod = "explicitRungeKutta4Delay";
    p.initialCondition  = "small";
    p.E_sp_det          = false;
    
    p.Ttrans = 1000.0;
    p.Teval  = 50.0;
    np.K_inj = 0.3;
    p.K_fb   = 0.25;
    p.phi    = 0.0;
    p.tau    = 8e-1;
    
    p.scanType = "downsweep";
    
    injectionBifurcationLineScan();
    
  }
  
  
  // bifurcation scan
  if (p.bifurcationScan == true) {
    
    p.networkType       = "singleLaser";
    p.timeSeriesMode    = "deterministic";
    p.integrationMethod = "explicitRungeKutta4";
    p.initialCondition  = "small";
    p.E_sp_det          = false;
    
    p.Ttrans = 150.0;
    p.Teval  = 50.0;
    
    // parameters for SNIPER bifurcation
    string bifurcationType = "SNIPER";
    np.K_inj               = 0.1;
    double DnuInjCrit      = -1.2588;
    
    // parameters for Hopf bifurcation
    //~ string bifurcationType = "Hopf";
    //~ np.K_inj               = 0.3;
    //~ double DnuInjCrit      = 2.647;
    
    bifurcationScan(bifurcationType, DnuInjCrit);
    
  }
  
  
  // feedback line-scan
  if (p.feedbackLineScan == true) {
    
    p.networkType       = "fourNodeAllToAll";
    p.timeSeriesMode    = "deterministic";
    p.integrationMethod = "explicitRungeKutta4Delay";
    p.initialCondition  = "synchron";
    p.E_sp_det          = false;
    p.Lyapunov          = false;
    
    p.Ttrans = 1000.0;
    p.Teval  = 50.0;
    p.tau    = 6e-1;
    
    p.scanType = "scan";
    
    feedbackLineScan();
    
  }
  
  
  // feedback bifurcation line-scan
  if (p.feedbackBifurcationLineScan == true) {
    
    p.networkType       = "singleLaserDelay";
    p.timeSeriesMode    = "deterministic";
    p.integrationMethod = "explicitRungeKutta4Delay";
    p.initialCondition  = "small";
    p.E_sp_det          = false;
    
    p.Ttrans = 1000.0;
    p.Teval  = 50.0;
    p.phi    = M_PI;
    p.tau    = 8e-1;
    
    p.scanType = "upsweep";
    
    feedbackBifurcationLineScan();
    
  }
  
  
  // rotating frame change
  if (p.rotatingFrameChange == true) {
    
    p.networkType       = "singleLaser";
    p.timeSeriesMode    = "deterministic";
    p.integrationMethod = "explicitRungeKutta4";
    p.initialCondition  = "small";
    p.E_sp_det          = true;
    
    p.Ttrans            = 0.0;
    p.Teval             = 10.0;
    
    rotatingFrameChange();
    
  }
  
  
  // find critical SNIPER bifurcation parameter DnuInj_crit
  if (p.findSNIPER == true) {
    
    // injection at K_inj=0.1:
    // beta = 5e-2: DnuInj = -1.2571471322458025
    // beta = 1e-2: DnuInj = -1.2584192267779066
    // beta = 5e-3: DnuInj = -1.2585501068308527
    // beta = 1e-3: DnuInj = -1.2586498268041225
    // beta = 5e-4: DnuInj = -1.2586619752501487
    // beta = 1e-4: DnuInj = -1.2586716429343596
    // beta = 1e-5: DnuInj = -1.2586738118937560
    
    // injection &   K_inj = 0.1, beta = 0       DnuInj = -1.8850417
    // feedback:     K_inj = 0.1, beta = 1e-3    DnuInj = -1.8850085
    //               K_inj = 0.1, beta = 1e-2    DnuInj = -1.8847083
    //               K_inj = 0.1, beta = 1e-1    DnuInj = -1.8815768
    //               K_inj = 0.1, beta = 1.0     DnuInj = -1.8372485
    
    // network at K_inj=0.1:
    // beta = 1.0:  DnuInj = -1.8372486054705173
    // beta = 5e-1: DnuInj = -1.8654959637641579
    // beta = 1e-1: DnuInj = -1.8815768550256819
    // beta = 5e-2: DnuInj = -1.8833448963022019
    // beta = 1e-2: DnuInj = -1.8847083704100058
    // beta = 5e-3: DnuInj = -1.8848754332408839
    // beta = 1e-3: DnuInj = -1.8850085277331479
    // beta = 5e-4: DnuInj = -1.8850251296055489
    // beta = 0:    DnuInj = -1.8850417236897399
    
    // network at K_inj=0.05:
    // beta = 1e-3: DnuInj = -1.2528700841245306
    
    
    TISI_Finder T;
    
    //~ p.networkType         = "singleLaser";
    //~ p.timeSeriesMode      = "deterministic";
    //~ p.integrationMethod   = "explicitRungeKutta4";
    //~ p.initialCondition    = "small";
    //~ p.E_sp_det            = true;
    //~ p.Ttrans              = 100.0;
    //~ p.Teval               = 250.0;
    //~ np.K_inj              = 0.1;
    
    //~ p.networkType         = "singleLaserDelay";
    //~ p.timeSeriesMode      = "deterministic";
    //~ p.integrationMethod   = "explicitRungeKutta4Delay";
    //~ p.initialCondition    = "small";
    //~ p.E_sp_det            = true;
    //~ p.Ttrans              = 1000.0;
    //~ p.Teval               = 250.0;
    //~ p.K_fb                = 0.05;
    //~ p.tau                 = 8e-2;
    //~ p.phi                 = M_PI / 2.0;
    //~ np.K_inj              = 0.1;
    
    p.networkType         = "fourNodeAllToAll";
    p.timeSeriesMode      = "deterministic";
    p.integrationMethod   = "explicitRungeKutta4Delay";
    p.initialCondition    = "synchron";
    p.E_sp_det            = true;
    p.Ttrans              = 1000.0;
    p.Teval               = 250.0;
    p.K_fb                = 0.05;
    p.tau                 = 8e-2;
    p.phi                 = M_PI / 2.0;
    np.K_inj              = 0.05;
    
    double bifurcationTol = 1e-15;
    double intensityTol   = 1e-2;
    bool   down           = true;
    
    T.findSNIPER(bifurcationTol, intensityTol, down);
    
  }
  
  
  // find critical homoclinic bifurcation parameter K_fb_crit
  if (p.findHomoclinic == true) {
    
    TISI_Finder T;
    
    p.networkType         = "singleLaserDelay";
    p.timeSeriesMode      = "deterministic";
    p.integrationMethod   = "explicitRungeKutta4Delay";
    p.initialCondition    = "small";
    p.E_sp_det            = true;
    
    p.Ttrans              = 1000.0;
    p.Teval               = 250.0;
    p.tau                 = 8e-1;
    p.phi                 = M_PI;
    
    double bifurcationTol = 1e-7;
    double intensityTol   = 1e-2;
    bool   left           = true;
    
    T.findHomoclinic(bifurcationTol, intensityTol, left);
    
  }
  
  
  // change of critical SNIPER bifurcation parameter DnuInj_crit with beta
  if (p.SNIPER_change == true) {
    
    p.networkType         = "singleLaserDelay";
    p.timeSeriesMode      = "deterministic";
    p.integrationMethod   = "explicitRungeKutta4Delay";
    p.initialCondition    = "small";
    p.E_sp_det            = true;
    
    p.Ttrans              = 1000.0;
    p.Teval               = 250.0;
    p.K_fb                = 0.05;
    p.tau                 = 8e-2;
    p.phi                 = M_PI / 2.0;
    np.K_inj              = 0.1;
    
    double bifurcationTol = 1e-15;
    double intensityTol   = 1e-2;
    bool   down           = true;
    
    SNIPER_change(bifurcationTol, intensityTol, down);
    
  }
  
  
  // change of critical homoclinic bifurcation parameter K_fb_crit with beta
  if (p.homoclinic_change == true) {
    
    p.networkType         = "singleLaserDelay";
    p.timeSeriesMode      = "deterministic";
    p.integrationMethod   = "explicitRungeKutta4Delay";
    p.initialCondition    = "small";
    p.E_sp_det            = true;
    
    p.Ttrans              = 1000.0;
    p.Teval               = 250.0;
    p.tau                 = 3e-1;
    p.phi                 = M_PI;
    
    double bifurcationTol = 1e-15;
    double intensityTol   = 1e-2;
    bool   left           = true;
    
    homoclinic_change(bifurcationTol, intensityTol, left);
    
  }
  
  
  // noise check
  if (p.noiseCheck == true) {
    
    p.networkType       = "singleLaser";
    p.dt                = 1e-4;
    int numberOfSamples = 1e6;
    noiseCheck(numberOfSamples);
    
  }
  
  
  // Tisi variance injection
  if (p.TisiVarianceInjection == true) {
    
    // single laser at K_inj=0.1
    //~ p.networkType         = "singleLaser";
    //~ p.timeSeriesMode      = "stochastic";
    //~ p.integrationMethod   = "explicitEulerMaruyama";
    //~ p.initialCondition    = "small";
    //~ p.E_sp_det            = false;
    //~ p.dt                  = 1e-4;
    //~ p.Ttrans              = 100.0;
    //~ p.Teval               = 250.0;
    //~ np.K_inj              = 0.1;
    
    
    // single laser with self feedback at K_inj=0.1, K_fb=0.05, phi=pi/2, tau=0.08
    //~ p.networkType         = "singleLaserDelay";
    //~ p.timeSeriesMode      = "stochastic";
    //~ p.integrationMethod   = "explicitEulerMaruyamaDelay";
    //~ p.initialCondition    = "small";
    //~ p.E_sp_det            = false;
    //~ p.dt                  = 1e-4;
    //~ p.Ttrans              = 1000.0;
    //~ p.Teval               = 250.0;
    //~ p.K_fb                = 0.05;
    //~ p.tau                 = 8e-2;
    //~ p.phi                 = M_PI / 2.0;
    //~ np.K_inj              = 0.1;
    
    // network at at K_inj=0.1, K_fb=0.05, phi=pi/2, tau=0.08
    p.networkType         = "fourNodeAllToAll";
    p.timeSeriesMode      = "stochastic";
    p.integrationMethod   = "explicitEulerMaruyamaDelay";
    p.initialCondition    = "synchron";
    p.E_sp_det            = false;
    p.dt                  = 1e-4;
    p.Ttrans              = 1000.0;
    p.Teval               = 250.0;
    p.K_fb                = 0.05;
    p.tau                 = 8e-2;
    p.phi                 = M_PI / 2.0;
    np.K_inj              = 0.1;
    
    int numberOfIntervals = 1000;
    
    TisiVarianceInjection(numberOfIntervals);
    
  }
  
  
  // Tisi variance feedback
  if (p.TisiVarianceFeedback == true) {
    
    p.networkType         = "singleLaserDelay";
    p.timeSeriesMode      = "stochastic";
    p.integrationMethod   = "explicitEulerMaruyamaDelay";
    p.initialCondition    = "small";
    p.E_sp_det            = false;
    
    p.Ttrans              = 1000.0;
    p.Teval               = 250.0;
    p.tau                 = 3e-2;
    p.phi                 = M_PI;
    
    int numberOfIntervals = 1000;
    
    TisiVarianceFeedback(numberOfIntervals);
    
  }
  
}

