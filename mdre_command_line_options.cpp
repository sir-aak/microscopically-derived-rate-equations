#include "mdre_parameters.hpp"
#include "mdre_command_line_options.hpp"


// enables command line options
void setCommandLineOptions (int argc, char* argv[])
{
  
  // link to instantiation of NodeParameters and Parameters class 
  // in mdre_processes.cpp
  extern NodeParameters np;
  extern Parameters p;
  
  if (argc > 0) {
    
    for (int i = 0; i < argc; i++) {
      
      // parameters for numerical integration
      
      if (string(argv[i]) == string("-dt") && argc > i+1) {
        p.dt         = atof(argv[i+1]);
        p.sqrt_dt    = sqrt(p.dt);
        p.Ito_dt     = sqrt(pow(p.dt, 3) / 3.0);
        p.transSteps = int(p.Ttrans / p.dt);
        p.evalSteps  = int(p.Teval / p.dt);
        p.nSteps     = int((p.Ttrans + p.Teval) / p.dt);
        p.n_tau      = int(p.tau / p.dt);
      }
      
      if (string(argv[i]) == string("-Ttrans") && argc > i+1) {
        p.Ttrans     = atof(argv[i+1]);
        p.transSteps = int(p.Ttrans / p.dt);
        p.nSteps     = int((p.Ttrans + p.Teval) / p.dt);
      }
      
      if (string(argv[i]) == string("-Teval") && argc > i+1) {
        p.Teval     = atof(argv[i+1]);
        p.evalSteps = int(p.Teval / p.dt);
        p.nSteps    = int((p.Ttrans + p.Teval) / p.dt);
      }
      
      if (string(argv[i]) == string("-outSteps") && argc > i+1) {
        p.outSteps = atoi(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-nSims") && argc > i+1) {
        p.nSims = atoi(argv[i+1]);
      }
      
      
      // method mode and type parameters
      
      if (string(argv[i]) == string("-networkType") && argc > i+1) {
        p.networkType = argv[i+1];
      }
      
      if (string(argv[i]) == string("-timeSeriesMode") && argc > i+1) {
        p.timeSeriesMode = argv[i+1];
      }
      
      if (string(argv[i]) == string("-integrationMethod") && argc > i+1) {
        p.integrationMethod = argv[i+1];
      }
      
      if (string(argv[i]) == string("-initialCondition") && argc > i+1) {
        p.initialCondition = argv[i+1];
      }
      
      if (string(argv[i]) == string("-scanType") && argc > i+1) {
        p.scanType = argv[i+1];
      }
      
      if (string(argv[i]) == string("-E_sp_det") && argc > i+1) {
        
        p.E_sp_det = atof(argv[i+1]);
        
        if (p.timeSeriesMode == "deterministic" && p.E_sp_det == false) {
          np.beta    = 0.0;
          cout << "\nbeta is set to zero. Set command line option -E_sp_det 1 for beta > 0.\n" << endl;
        }
        
      }
      
      if (string(argv[i]) == string("-Lyapunov") && argc > i+1) {
        p.Lyapunov = atof(argv[i+1]);
      }
      
      
      // control parameters:
      
      if (string(argv[i]) == string("-J") && argc > i+1) {
        np.J = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-sc") && argc > i+1) {
        np.sc = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-beta") && argc > i+1) {
        np.beta = atof(argv[i+1]);
      }
      
      
      // injection parameters:
      
      if (string(argv[i]) == string("-K_inj") && argc > i+1) {
        np.K_inj = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-DnuInj") && argc > i+1) {
        np.DnuInj = atof(argv[i+1]);
      }
      
      
      // feedback coupling parameters
      
      if (string(argv[i]) == string("-K_fb") && argc > i+1) {
        p.K_fb = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-phi") && argc > i+1) {
        p.phi = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-tau") && argc > i+1) {
        p.tau   = atof(argv[i+1]);
        p.n_tau = int(p.tau / p.dt);
      }
      
      
      // laser parameters:
      
      if (string(argv[i]) == string("-a_L") && argc > i+1) {
        np.a_L = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-f_act") && argc > i+1) {
        np.f_act   = atof(argv[i+1]);
        np.f_inact = 1.0 - np.f_act;
      }
      
      if (string(argv[i]) == string("-g_GS") && argc > i+1) {
        np.g_GS = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-N_QD") && argc > i+1) {
        np.N_QD = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-R_w_loss") && argc > i+1) {
        np.R_w_loss = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-T_e_eq") && argc > i+1) {
        np.T_e_eq             = atof(argv[i+1]);
        np.factor_S_e_rel_out = exp((np.eps_GS_e - np.eps_ES_e) 
                              / (np.kB * np.T_e_eq));
      }
      
      if (string(argv[i]) == string("-T_h_eq") && argc > i+1) {
        np.T_h_eq             = atof(argv[i+1]);
        np.factor_S_h_rel_out = exp((np.eps_GS_h - np.eps_ES_h) 
                              / (np.kB * np.T_h_eq));
      }
      
      if (string(argv[i]) == string("-W_G") && argc > i+1) {
        np.W_GS = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-W_ES") && argc > i+1) {
        np.W_ES = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-kappa") && argc > i+1) {
        np.kappa = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-delta_omega_ES") && argc > i+1) {
        np.delta_omega_ES = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-delta_omega_e_QW") && argc > i+1) {
        np.delta_omega_e_QW = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-delta_omega_h_QW") && argc > i+1) {
        np.delta_omega_h_QW = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-eps_ES_e") && argc > i+1) {
        np.eps_ES_e           = atof(argv[i+1]);
        np.factor_S_e_rel_out = exp((np.eps_GS_e - np.eps_ES_e) 
                              / (np.kB * np.T_e_eq));
      }
      
      if (string(argv[i]) == string("-eps_GS_e") && argc > i+1) {
        np.eps_GS_e           = atof(argv[i+1]);
        np.factor_S_e_rel_out = exp((np.eps_GS_e - np.eps_ES_e) 
                              / (np.kB * np.T_e_eq));
      }
      
      if (string(argv[i]) == string("-eps_GS_h") && argc > i+1) {
        np.eps_GS_h           = atof(argv[i+1]);
        np.factor_S_h_rel_out = exp((np.eps_GS_h - np.eps_ES_h) 
                              / (np.kB * np.T_h_eq));
      }
      
      if (string(argv[i]) == string("-eps_ES_h") && argc > i+1) {
        np.eps_ES_h           = atof(argv[i+1]);
        np.factor_S_h_rel_out = exp((np.eps_GS_h - np.eps_ES_h) 
                              / (np.kB * np.T_h_eq));
      }
      
      
      // fit parameters for charge-carrier scattering processes:
      
      if (string(argv[i]) == string("-A_GS_e") && argc > i+1) {
        np.A_GS_e = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-B_GS_e") && argc > i+1) {
        np.B_GS_e = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-A_GS_h") && argc > i+1) {
        np.A_GS_h = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-B_GS_h") && argc > i+1) {
        np.B_GS_h = atof(argv[i+1]);
      }
      
      
      if (string(argv[i]) == string("-A_ES_e") && argc > i+1) {
        np.A_ES_e = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-B_ES_e") && argc > i+1) {
        np.B_ES_e = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-A_ES_h") && argc > i+1) {
        np.A_ES_h = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-B_ES_h") && argc > i+1) {
        np.B_ES_h = atof(argv[i+1]);
      }
      
      
      if (string(argv[i]) == string("-C_e") && argc > i+1) {
        np.C_e = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-D_e") && argc > i+1) {
        np.D_e = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-C_h") && argc > i+1) {
        np.C_h = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-D_h") && argc > i+1) {
        np.D_h = atof(argv[i+1]);
      }
      
      
      // function switches
      
      if (string(argv[i]) == string("-timeSeries") && argc > i+1) {
        p.timeSeries = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-inputOutput") && argc > i+1) {
        p.inputOutput = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-nonlinearScatteringRates") 
          && argc > i+1) {
        p.nonlinearScatteringRates = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-scatteringProcesses") 
          && argc > i+1) {
        p.scatteringProcesses = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-injectionLineScan") && argc > i+1) {
        p.injectionLineScan = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-injectionBifurcationLineScan") && argc > i+1) {
        p.injectionBifurcationLineScan = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-bifurcationScan") && argc > i+1) {
        p.bifurcationScan = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-feedbackLineScan") && argc > i+1) {
        p.feedbackLineScan = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-feedbackBifurcationLineScan") 
          && argc > i+1) {
        p.feedbackBifurcationLineScan = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-rotatingFrameChange") && argc > i+1) {
        p.rotatingFrameChange = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-findSNIPER") && argc > i+1) {
        p.findSNIPER = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-findHomoclinic") && argc > i+1) {
        p.findHomoclinic = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-SNIPER_change") && argc > i+1) {
        p.SNIPER_change = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-homoclinic_change") && argc > i+1) {
        p.homoclinic_change = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-noiseCheck") && argc > i+1) {
        p.noiseCheck = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-TisiVarianceInjection") && argc > i+1) {
        p.TisiVarianceInjection = atof(argv[i+1]);
      }
      
      if (string(argv[i]) == string("-TisiVarianceFeedback") && argc > i+1) {
        p.TisiVarianceFeedback = atof(argv[i+1]);
      }
      
    }
    
  }
  
}

