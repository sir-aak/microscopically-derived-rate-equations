#include "mdre_parameters.hpp"
#include "mdre_node_equations.hpp"
#include "mdre_network.hpp"
#include "mdre_solver.hpp"
#include "mdre_time_series.hpp"
#include "mdre_auxiliary.hpp"


// link to Parameters instantiation in mdre_processes.cpp
extern NodeParameters np;
extern Parameters p;


// adjusts numeric parameters to new dt, Ttrans, Teval and tau
void updateNumericParameters ()
{
  
  p.sqrt_dt       = sqrt(p.dt / 2.0);
  p.nSteps        = int((p.Ttrans + p.Teval) / p.dt);
  p.transSteps    = int(p.Ttrans / p.dt);
  p.evalSteps     = int(p.Teval / p.dt);
  p.n_tau         = int(p.tau / p.dt);
  
}



// returns vector for storage of nonlinear scattering rates
arma::rowvec nonlinearScatteringRatesVector (TimeSeries& TS, double w)
{
  
  TS.N.nodes[0].w_e = w;
  TS.N.nodes[0].w_h = w;
  
  arma::rowvec nonlinearScatteringData = {
    w,
    TS.N.nodes[0].S_GS_e_cap_in(),
    TS.N.nodes[0].S_GS_h_cap_in(),
    TS.N.nodes[0].S_ES_e_cap_in(),
    TS.N.nodes[0].S_ES_h_cap_in(),
    TS.N.nodes[0].S_e_rel_in(),
    TS.N.nodes[0].S_h_rel_in(),
    
    TS.N.nodes[0].S_GS_e_cap_out(),
    TS.N.nodes[0].S_GS_h_cap_out(),
    TS.N.nodes[0].S_ES_e_cap_out(),
    TS.N.nodes[0].S_ES_h_cap_out(),
    TS.N.nodes[0].S_e_rel_out(),
    TS.N.nodes[0].S_h_rel_out()
  };
  
  return (nonlinearScatteringData);
  
}



// returns vector for storage of scattering processes
arma::rowvec scatteringProcessesVector (TimeSeries& TS, double w)
{
  
  TS.N.nodes[0].w_e = w;
  TS.N.nodes[0].w_h = w;
  
  arma::rowvec scatteringProcessesData = {
    w,
    TS.N.nodes[0].S_GS_e_cap_act(),
    TS.N.nodes[0].S_GS_h_cap_act(),
    TS.N.nodes[0].S_GS_e_cap_inact(),
    TS.N.nodes[0].S_GS_h_cap_inact(),
    TS.N.nodes[0].S_ES_e_cap(),
    TS.N.nodes[0].S_ES_h_cap(),
    TS.N.nodes[0].S_e_rel_act(),
    TS.N.nodes[0].S_h_rel_act(),
    TS.N.nodes[0].S_e_rel_inact(),
    TS.N.nodes[0].S_h_rel_inact()
  };
  
  return (scatteringProcessesData);
  
}



// returns inverse electron scattering lifetime
double tau_e_inverse (double rho_ES_e, double w_e)
{
  
  updateNumericParameters();
  TimeSeries TS;
  
  return (1.0 / (TS.N.nodes[0].S_GS_e_cap_in() 
  + TS.N.nodes[0].S_GS_e_cap_out() 
  + TS.N.nodes[0].S_e_rel_in() * TS.N.nodes[0].rho_ES_e 
  + TS.N.nodes[0].S_e_rel_out() * (1.0 - TS.N.nodes[0].rho_ES_e)));
  
}



// returns inverse hole scattering lifetime
double tau_h_inverse (double rho_ES_h, double w_h)
{
  
  updateNumericParameters();
  TimeSeries TS;
  
  return (1.0 / (TS.N.nodes[0].S_GS_h_cap_in() 
  + TS.N.nodes[0].S_GS_h_cap_out() 
  + TS.N.nodes[0].S_h_rel_in() * TS.N.nodes[0].rho_ES_h 
  + TS.N.nodes[0].S_h_rel_out() * (1.0 - TS.N.nodes[0].rho_ES_h)));
  
}



// creates argument file for ..
void createArgumentFile (string argfileType)
{
  
  ofstream file;
  
  // .. 2D injection bifurcation diagram
  // needed when using injectionLineScan() on condor
  if (argfileType == "injection") {
    
    double K_injStart = 0.0;
    double K_injEnd   = 0.5;
    double K_injStep  = pow(2.0, -9);
    
    file.open("injection_argfile.txt");
    
    for (double K_inj = K_injStart; K_inj <= K_injEnd; K_inj += K_injStep) {
      file << "-injectionLineScan 1 " << "-K_inj " << K_inj << endl;
    }
    
  }
  
  // .. K_fb-phi feedback plot
  // needed when using feedbackLineScan() on condor
  else if (argfileType == "feedback") {
    
    double K_fbStart     = 0.0;
    double K_fbEnd       = 0.5;
    double K_fbIncrement = pow(2.0, -9);
    
    file.open("feedback_argfile.txt");
    
    for (double K_fb = K_fbStart; K_fb <= K_fbEnd; 
         K_fb += K_fbIncrement) {
      
      file << "-feedbackLineScan 1 -K_fb " << K_fb << endl;
      
    }
    
  }
  
  // .. DnuInj_crit(beta) plot
  // needed when using SNIPER_change() on condor
  else if (argfileType == "SNIPER_change") {
    
    double betaExpStart     = -7.0;
    double betaExpEnd       = 0.0;
    double betaExpIncrement = 0.0625;
    
    file.open("SNIPER_argfile.txt");
    
    for (double betaExp = betaExpStart; betaExp <= betaExpEnd; 
         betaExp += betaExpIncrement) {
      
      file << "-SNIPER_change 1 " 
           << "-beta " << pow(10.0, betaExp) << endl;
      
    }
    
  }
  
  // .. K_fb_crit(beta) plot
  // needed when using homoclinic_change() on condor
  else if (argfileType == "homoclinic_change") {
    
    double betaExpStart     = -7.0;
    double betaExpEnd       = 0.0;
    double betaExpIncrement = 0.0625;
    
    file.open("homoclinic_argfile.txt");
    
    for (double betaExp = betaExpStart; betaExp <= betaExpEnd; 
         betaExp += betaExpIncrement) {
      
      file << "-homoclinic_change 1 " 
           << "-beta " << pow(10.0, betaExp) << endl;
      
    }
    
  }
  
  // .. coherence resonance diagram for injection
  // needed when using TisiVarianceInjection() on condor
  else if (argfileType == "coherenceResonanceInjection") {
    
    double betaExpStart     = -4.0; //-3.375;
    double betaExpEnd       = 0.0;
    double betaExpIncrement = 0.03125;
    
    file.open("coherence_resonance_injection_argfile.txt");
    
    for (double betaExp = betaExpStart; betaExp <= betaExpEnd; 
         betaExp += betaExpIncrement) {
      
      file << "-TisiVarianceInjection 1 " 
           << "-beta " << pow(10.0, betaExp) << endl;
      
    }
    
  }
  
  // .. coherence resonance diagram for injection with delayed feedback
  // needed when using TisiVarianceInjection() on condor
  else if (argfileType == "coherenceResonanceInjectionFeedback") {
    
    double betaExpStart     = -6.0; // -2.4375;
    double betaExpEnd       = 0.0;
    double betaExpIncrement = 0.03125;
    
    file.open("coherence_resonance_injection_feedback_argfile.txt");
    
    for (double betaExp = betaExpStart; betaExp <= betaExpEnd; 
         betaExp += betaExpIncrement) {
      
      file << "-TisiVarianceInjection 1 " 
           << "-beta " << pow(10.0, betaExp) << endl;
      
    }
    
  }
  
  // .. coherence resonance diagram for network
  // needed when using TisiVarianceInjection() on condor
  else if (argfileType == "coherenceResonanceNetwork") {
    
    double betaExpStart     = 0.0; //-5.5;
    double betaExpEnd       = -12.0;
    double betaExpIncrement = -0.03125;
    
    file.open("coherence_resonance_network_argfile.txt");
    
    for (double betaExp = betaExpStart; betaExp >= betaExpEnd; 
         betaExp += betaExpIncrement) {
      
      file << "-TisiVarianceInjection 1 " 
           << "-beta " << pow(10.0, betaExp) << endl;
      
    }
    
  }
  
  else {
    cerr << "choose argument file type in {injection, feedback, SNIPER_change, homoclinic_change, coherenceResonanceInjection, coherenceResonanceInjectionFeedback, coherenceResonanceFeedback}" << endl;
    exit(0);
  }
  
  file.close();
  
}

