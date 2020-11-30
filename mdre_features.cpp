#include "mdre_parameters.hpp"
#include "mdre_node_equations.hpp"
#include "mdre_network.hpp"
#include "mdre_solver.hpp"
#include "mdre_time_series.hpp"
#include "mdre_auxiliary.hpp"
#include "mdre_analyzer.hpp"
#include "mdre_TISI_finder.hpp"
#include "mdre_features.hpp"


// link to instantiation of NodeParameters and Parameters class 
// in mdre_processes.cpp
extern NodeParameters np;
extern Parameters p;


// generates data for deterministic or stochastic time series
void timeSeries ()
{
  
  // start of time measure
  auto t_start = chrono::steady_clock::now();
  
  updateNumericParameters();
  TimeSeries TS;
  
  // define path and file name for time series file
  
  string path;
  
  if (p.timeSeriesMode == "deterministic") {
    path = "time_series/time_series.txt";
  }
  
  else if (p.timeSeriesMode == "stochastic") {
    path = "time_series/time_series_stochastic.txt";
  }
  
  ofstream file(path);
  file.precision(12);
  
  // simulations loop
  for (int i = 0; i < p.nSims; i++) {
    
    TS.t = 0.0;
    
    // set initial condition
    TS.N.setInitialConditions();
    
    // write initial condition in time series file
    if (p.Ttrans == 0.0) TS.writeTimeSeries(file);
    
    // time series loop
    for (long long int j = 1; j < p.nSteps; j++) {
      
      // time series step
      TS.timeSeriesStep(TS.N, j);
      
      // waiting for the transient to evaluate
      if (j >= p.transSteps && j % p.outSteps == 0) {
        TS.writeTimeSeries(file);
      }
      
      // pertubation on electric fields of lasers in network
      // to check their synchronization
      if (TS.t >= 50.0 && TS.t < 55.0 
          && p.networkType == "fourNodeAllToAll" 
          && p.initialCondition == "synchron") {
        TS.networkPerturbation();
      }
      
      // progress indicator
      if (p.nSims == 1 && j % int(p.nSteps / 100.0) == 0) {
        cout << int(j * 100.0 / p.nSteps) << " %" << endl;
      }
      
    }
    
    cout << "simulation " << i+1 << " of " << p.nSims << endl;
    
  }
  
  file.close();
  
  // end of time measure
  auto t_end = chrono::steady_clock::now();
  
  if (p.timeSeriesMode == "deterministic") {
    cout << "\nDeterministic time series was computed with the following parameters:\n" << endl;
  }
  
  else if (p.timeSeriesMode == "stochastic") {
    cout << "\nStochastic time series was computed with the following parameters:\n" << endl;
  }
  
  TS.showParameters();
  
  cout << "\ncomputation time: " 
       << chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count() / 1e3 
       << " s\n" << endl;
  
}



// generates data for input-output bifurcation diagram
void inputOutput_mdre ()
{
  
  // start of time measure
  auto t_start = chrono::steady_clock::now();
  
  updateNumericParameters();
  TimeSeries TS;
  
  // define path and file name for input-output file
  ofstream file("input_output/input_output_mdre.txt");
  file.precision(12);
  
  double JStart = 0.0;
  double JEnd   = 50.0;
  double JStep  = 0.0625;
  
  int JSteps = int((JEnd - JStart) / JStep);
  
  cout << JSteps << endl;
  
  np.J = JStart;
  
  // pump current loop
  for (int j = 0; j <= JSteps; j++) {
    
    TS.N.setInitialConditions();
    TS.t = 0.0;
    
    // time series loop
    for (long long int i = 1; i < p.nSteps; i++) {
      TS.timeSeriesStep(TS.N, i);
    }
    
    file << np.J                         << '\t' 
         << norm(TS.N.nodes[0].E)        << '\t' 
         << TS.N.nodes[0].rho_GS_e_act   << '\t' 
         << TS.N.nodes[0].rho_GS_h_act   << '\t' 
         << TS.N.nodes[0].rho_GS_e_inact << '\t' 
         << TS.N.nodes[0].rho_GS_h_inact << '\t' 
         << TS.N.nodes[0].rho_ES_e       << '\t' 
         << TS.N.nodes[0].rho_ES_h       << '\t' 
         << TS.N.nodes[0].w_e            << '\t' 
         << TS.N.nodes[0].w_h            << endl;
    
    cout << "J = " << np.J << '\t' << "rho_GS_e_act + rho_GS_h_act = " 
         << TS.N.nodes[0].rho_GS_e_act + TS.N.nodes[0].rho_GS_h_act << endl;
    
    np.J += JStep;
    
  }
  
  file.close();
  
  // end of time measure
  auto t_end = chrono::steady_clock::now();
  
  cout << "\nInput-output was computed with the following parameters:\n" << endl;
  
  TS.showParameters();
  
  cout << "\ncomputation time: " 
       << chrono::duration_cast<chrono::seconds>(t_end - t_start).count() / 60.0
       << " min\n" << endl;
  
}



// generates data for nonlinear scattering rates
void nonlinearScatteringRates ()
{
  
  updateNumericParameters();
  TimeSeries TS;
  
  double wStart = 0.0;
  double wEnd   = 10.0;
  double wStep  = 0.0625;
  
  double w  = wStart;
  int nRows = int((wEnd - wStart) / wStep);
  int nCols = nonlinearScatteringRatesVector(TS, wStart).n_elem;
  
  arma::mat data(nRows, nCols);
  data.zeros();
  
  for (int i = 0; i < nRows; i++) {
    data.row(i) = nonlinearScatteringRatesVector(TS, w);
    w          += wStep;
  }
  
  // define path and file name for nonlinear scattering rates file
  ofstream file("scattering_processes/nonlinear_scattering_rates_mdre.txt");
  file.precision(12);
  file << data << endl;
  file.close();
  
}



// generates data for scattering processes between charge-carrier states
void scatteringProcesses ()
{
  
  updateNumericParameters();
  TimeSeries TS;
  
  double wStart = 0.0;
  double wEnd   = 10.0;
  double wStep  = 0.0625;
  
  double  w = wStart;
  int nRows = int((wEnd - wStart) / wStep);
  int nCols = scatteringProcessesVector(TS, wStart).n_elem;
  
  arma::mat data(nRows, nCols);
  data.zeros();
  
  p.Ttrans = 0.0;
  p.Teval  = 10.0;
  p.nSteps = int((p.Ttrans + p.Teval) / p.dt);
  
  for (int i = 1; i < p.nSteps; i++) {
    TS.timeSeriesStep(TS.N, i);
  }
  
  for (int i = 0; i < nRows; i++) {
    data.row(i) = scatteringProcessesVector(TS, w);
    w          += wStep;
  }
  
  // define path and file name for scattering processes file
  ofstream file("scattering_processes/scattering_processes_mdre.txt");
  file.precision(12);
  file << data << endl;
  file.close();
  
}



// generates data for number of unique maxima and period diagram
// in single laser case also Lyapunov exponent diagram
// in network case also network state diagram
// along DnuInj at given injection strength K_inj
// use argfileType "injection"
void injectionLineScan ()
{
  
  // start of time measure
  auto t_start = chrono::steady_clock::now();
  
  updateNumericParameters();
  TimeSeries TS;
  Analyzer   A;
  
  double DnuInj;
  double DnuInjStart;
  double DnuInjEnd;
  double DnuInjStep;
  
  // define path and file name for injection line scan file
  string path = "injection/1/injection_line_scan_mdre_K_inj=" 
              + to_string(np.K_inj) + ".txt";
  ofstream file(path); file.precision(12);
  
  if (p.scanType == "scan" || p.scanType == "upsweep") {
    DnuInjStart = -6.0;
    DnuInjEnd   =  6.0;
    DnuInjStep  =  0.03125;
  }
  
  else if (p.scanType == "downsweep") {
    DnuInjStart =  6.0;
    DnuInjEnd   = -6.0;
    DnuInjStep  = -0.03125;
  }
  
  else {
    cerr << "\nChoose scan type from {scan, upsweep, downsweep}.\n" << endl;
    exit(0);
  }
  
  cout << "K_inj = " << np.K_inj << endl;
  
  // detuning loop, <= for upsweep, >= for downsweep
  for (DnuInj = DnuInjStart; DnuInj <= DnuInjEnd; DnuInj += DnuInjStep) {
    
    TS.t      = 0.0;
    np.DnuInj = DnuInj;
    
    // set initial state of system for new detuning parameter
    if (p.scanType == "scan") TS.N.setInitialConditions();
    
    // time series loop
    for (long long int j = 1; j < p.nSteps; j++) {
      
      // perform time series step
      TS.timeSeriesStep(TS.N, j);
      
      // pertubation to check if solution stays synchronous
      if (TS.t >= 50.0 && TS.t < 55.0 
          && p.networkType == "fourNodeAllToAll" 
          && p.initialCondition == "synchron") {
        TS.networkPerturbation();
      }
      
      // store time, intensity and Lyapunov values after transient time
      if (j >= p.transSteps) {
        
        TS.timeValues(j - p.transSteps)       = TS.t;
        TS.intensityValues0(j - p.transSteps) = norm(TS.N.nodes[0].E);
        
        if ((p.networkType == "singleLaser" 
          || p.networkType == "singleLaserDelay") && p.Lyapunov == true) {
          TS.LyapunovValues(j - p.transSteps) = TS.N.LyapunovExponent();
        }
        
        if (p.networkType == "fourNodeAllToAll") {
          TS.intensityValues1(j - p.transSteps) = norm(TS.N.nodes[1].E);
          TS.intensityValues2(j - p.transSteps) = norm(TS.N.nodes[2].E);
          TS.intensityValues3(j - p.transSteps) = norm(TS.N.nodes[3].E);
        }
        
      }
      
    }
    
    // analyzing time series data
    
    if (p.networkType == "singleLaser" 
     || p.networkType == "singleLaserDelay") {
      
      A.fixedPointCheck(TS.intensityValues0);
      
      if (p.Lyapunov == true) {
        A.analyze(TS.timeValues, TS.intensityValues0, TS.LyapunovValues);
      }
      
      else A.analyze(TS.timeValues, TS.intensityValues0);
      
      A.writeInjectionSingleLaserData(file);
      
    }
    
    else if (p.networkType == "fourNodeAllToAll") {
      A.fixedPointCheck(TS.intensityValues0, TS.intensityValues1, 
                        TS.intensityValues2, TS.intensityValues3);
      A.analyze(TS.timeValues, TS.intensityValues0, TS.intensityValues1, 
                TS.intensityValues2, TS.intensityValues3);
      A.writeInjectionLaserNetworkData(file);
    }
    
    A.resetData();
    
  }
  
  file.close();
  
  // end of time measure
  auto t_end = chrono::steady_clock::now();
  
  cout << "\nInjection line scan was computed with the following parameters:\n" << endl;
  
  TS.showParameters();
  
  cout << "\ncomputation time: " 
       << chrono::duration_cast<chrono::seconds>(t_end - t_start).count() / 60.0
       << " min\n" << endl;
  
}



// generates data for bifurcation line scan along master-slave detuning
// DnuInj at given injection coupling strength K_inj
// has to be called with singleLaser or singleLaserDelay mode and with 
// command line option -K_inj
void injectionBifurcationLineScan ()
{
  
  // start of time measure
  auto t_start = chrono::steady_clock::now();
  
  updateNumericParameters();
  TimeSeries TS;
  Analyzer A;
  
  double DnuInj;
  double DnuInjStart;
  double DnuInjEnd;
  double DnuInjStep;
  
  // define path and file name for feedback line scan file
  string path = "injection/bifurcation_line_scan/injection_bifurcation_line_scan_K_inj="
              + to_string(np.K_inj) + "_K_fb=" + to_string(p.K_fb) 
              + "_phi=" + to_string(p.phi) + "_tau=" + to_string(p.tau);
  if (p.scanType == "upsweep") path += "_upsweep";
  else if (p.scanType == "downsweep") path += "_downsweep";
  path += ".txt";
  ofstream file(path); file.precision(12);
  
  if (p.scanType == "scan" || p.scanType == "upsweep") {
    DnuInjStart = -6.0;
    DnuInjEnd   =  6.0;
    DnuInjStep  =  0.03125;
  }
  
  else if (p.scanType == "downsweep") {
    DnuInjStart =  6.0;
    DnuInjEnd   = -6.0;
    DnuInjStep  = -0.03125;
  }
  
  else {
    cerr << "\nChoose scan type from {scan, upsweep, downsweep}.\n" << endl;
    exit(0);
  }
  
  cout << "K_inj = " << np.K_inj << endl;
  
  // detuning loop, <= for upsweep, >= for downsweep
  for (DnuInj = DnuInjStart; DnuInj >= DnuInjEnd; DnuInj += DnuInjStep) {
    
    np.DnuInj = DnuInj;
    
    TS.t = 0.0;
    
    // set initial state of system for new detuning parameter
    if (p.scanType == "scan") TS.N.setInitialConditions();
    
    // loop for creation of time series data
    for (long long int j = 1; j < p.nSteps; j++) {
      
      // perform time series step
      TS.timeSeriesStep(TS.N, j);
      
      // store time and intensity values
      if (j >= p.transSteps) {
        TS.timeValues(j - p.transSteps)       = TS.t;
        TS.intensityValues0(j - p.transSteps) = norm(TS.N.nodes[0].E);
      }
      
    }
    
    A.fixedPointCheck(TS.intensityValues0);
    A.analyze(TS.timeValues, TS.intensityValues0);
    A.writeInjectionBifurcationData(file);
    A.resetData();
    
  }
  
  file.close();
  
  // end of time measure
  auto t_end = chrono::steady_clock::now();
  
  cout << "\nInjection bifurcation line scan was computed with the following parameters:\n" << endl;
  
  TS.showParameters();
  
  cout << "\ncomputation time: " 
       << chrono::duration_cast<chrono::seconds>(t_end - t_start).count() / 60.0
       << " min\n" << endl;
  
}



// generates data for SNIPER or Hopf bifurcation scan
void bifurcationScan (string bifurcationType, double DnuInjCrit)
{
  
  // start of time measure
  auto t_start = chrono::steady_clock::now();
  
  updateNumericParameters();
  TimeSeries TS;
  Analyzer A;
  
  double mu;				// mu = DnuInjcrit - DnuInj
  double greatestMaximum;
  
  // define path and file name for bifurcation file
  
  string path;
  
  if (bifurcationType == "SNIPER")    path = "bifurcation_scans/sniper.txt";
  else if (bifurcationType == "Hopf") path = "bifurcation_scans/hopf.txt";
  else {
    cerr << "\nChoose bifurcation type from {SNIPER, Hopf}.\n" << endl;
    exit(0);
  }
  
  ofstream file(path); file.precision(12);
  
  // loop for bifurcation parameter mu
  for (double muExp = -3.0; muExp <= -1.0; muExp += 0.125) {
    
    mu = pow(10.0, muExp);
    
    if (bifurcationType == "SNIPER")    np.DnuInj = DnuInjCrit - mu;
    else if (bifurcationType == "Hopf") np.DnuInj = DnuInjCrit + mu;
    
    // set initial state of system
    TS.t = 0.0;
    TS.N.setInitialConditions();
    
    // loop for creation of time series data
    for (long long int i = 1; i < p.nSteps; i++) {
      
      // perform time series step
      TS.timeSeriesStep(TS.N, i);
      
      // store time and intensity values
      if (i >= p.transSteps) {
        TS.timeValues(i - p.transSteps)       = TS.t;
        TS.intensityValues0(i - p.transSteps) = norm(TS.N.nodes[0].E);
      }
      
    }
    
    // analyzing time series data
    A.analyze(TS.timeValues, TS.intensityValues0);
    
    // find greatest maximum and smallest minimum
    greatestMaximum = *max_element(A.maxima0.begin(), A.maxima0.end());
    
    file << mu << '\t' << A.period0 << '\t' << greatestMaximum << endl;
    
    cout << "mu = "        << mu              << '\t' 
         << "period = "    << A.period0       << '\t' 
         << "amplitude = " << greatestMaximum << '\t' 
         << "muExp = "     << muExp           << endl;
    
    A.resetData();
    
  }
  
  file.close();
  
  // end of time measure
  auto t_end = chrono::steady_clock::now();
  
  cout << "\nbifurcation scan was computed with the following parameters:" << endl;
  
  TS.showParameters();
  
  cout << "\ncomputation time: " 
       << chrono::duration_cast<chrono::seconds>(t_end - t_start).count() 
       << " s\n" << endl;
  
}



// generates data for K_fb-phi diagram along feedback coupling phase phi 
// at given feedback coupling strength K_fb
// has to be called with singleLaserDelay or fourNodeAllToAll mode 
// and with command line option -K_fb
// use argfileType "feedback"
void feedbackLineScan ()
{
  
  // start of time measure
  auto t_start = chrono::steady_clock::now();
  
  updateNumericParameters();
  TimeSeries TS;
  Analyzer A;
  
  int i;
  int iStart;
  int iEnd;
  
  // define path and file name for feedback line scan file
  string path = "feedback/feedback_line_scan_mdre_K_fb=" 
              + to_string(p.K_fb) + ".txt";
  ofstream file(path); file.precision(12);
  
  if (p.scanType == "scan" || p.scanType == "upsweep") {
    iStart   = 0;
    iEnd     = 256;
  }
  
  else if (p.scanType == "downsweep") {
    iStart   = 256;
    iEnd     = 0;
  }
  
  else {
    cerr << "\nChoose scan type from {scan, upsweep, downsweep}.\n" << endl;
    exit(0);
  }
  
  cout << "K_fb = " << p.K_fb << endl;
  
  // coupling phase loop, <= for upsweep, >= for downsweep
  for (i = iStart; i <= iEnd; ) {
    
    p.phi = i * 2.0 * M_PI / 256.0;
    
    if (p.scanType == "scan" || p.scanType == "upsweep") i++;
    else if (p.scanType == "downsweep") i--;
    
    TS.t = 0.0;
    
    // set initial state of system for new detuning parameter
    if (p.scanType == "scan") TS.N.setInitialConditions();
    
    // time series loop
    for (long long int j = 1; j < p.nSteps; j++) {
      
      // perform time series step
      TS.timeSeriesStep(TS.N, j);
      
      // pertubation to check if solution stays synchronous
      if (TS.t >= 50.0 && TS.t < 55.0 
          && p.networkType == "fourNodeAllToAll" 
          && p.initialCondition == "synchron") {
        TS.networkPerturbation();
      }
      
      // store time and intensity values
      if (j >= p.transSteps) {
        
        TS.timeValues(j - p.transSteps)       = TS.t;
        TS.intensityValues0(j - p.transSteps) = norm(TS.N.nodes[0].E);
        
        if (p.networkType == "singleLaserDelay" && p.Lyapunov == true) {
          TS.LyapunovValues(j - p.transSteps) = TS.N.LyapunovExponent();
        }
        
        if (p.networkType == "fourNodeAllToAll") {
          TS.intensityValues1(j - p.transSteps) = norm(TS.N.nodes[1].E);
          TS.intensityValues2(j - p.transSteps) = norm(TS.N.nodes[2].E);
          TS.intensityValues3(j - p.transSteps) = norm(TS.N.nodes[3].E);
        }
        
      }
      
    }
    
    // analyzing time series data
    
    if (p.networkType == "singleLaserDelay") {
      
      A.fixedPointCheck(TS.intensityValues0);
      
      if (p.Lyapunov == true) {
        A.analyze(TS.timeValues, TS.intensityValues0, TS.LyapunovValues);
      }
      
      else A.analyze(TS.timeValues, TS.intensityValues0);
      
      A.writeFeedbackSingleLaserData(file);
      
    }
    
    else if (p.networkType == "fourNodeAllToAll") {
      A.fixedPointCheck(TS.intensityValues0, TS.intensityValues1, 
                        TS.intensityValues2, TS.intensityValues3);
      A.analyze(TS.timeValues, TS.intensityValues0, TS.intensityValues1, 
                TS.intensityValues2, TS.intensityValues3);
      A.writeFeedbackLaserNetworkData(file);
    }
    
    A.resetData();
    
  }
  
  file.close();
  
  // end of time measure
  auto t_end = chrono::steady_clock::now();
  
  cout << "\nFeedback line scan was computed with the following parameters:\n" << endl;
  
  TS.showParameters();
  
  cout << "\ncomputation time: " 
       << chrono::duration_cast<chrono::seconds>(t_end - t_start).count() / 60.0
       << " min\n" << endl;
  
}



// generates data for bifurcation line scan along feedback coupling
// strength K_fb at given feedback coupling phase phi
// has to be called with singleLaser or singleLaserDelay mode and with 
// command line option -phi
void feedbackBifurcationLineScan ()
{
  
  // start of time measure
  auto t_start = chrono::steady_clock::now();
  
  updateNumericParameters();
  TimeSeries TS;
  Analyzer A;
  
  double K_fb;
  double K_fbStart;
  double K_fbEnd;
  double K_fbStep;
  
  // define path and file name for feedback line scan file
  string path = "feedback/bifurcation_line_scan/feedback_bifurcation_line_scan_phi="
              + to_string(p.phi) + "_tau=" + to_string(p.tau);
  if (p.scanType == "upsweep") path += "_upsweep";
  else if (p.scanType == "downsweep") path += "_downsweep";
  path += ".txt";
  ofstream file(path); file.precision(12);
  
  if (p.scanType == "scan" || p.scanType == "upsweep") {
    K_fbStart = 0.0;
    K_fbEnd   = 0.5;
    K_fbStep  = 0.0009765625;
  }
  
  else if (p.scanType == "downsweep") {
    K_fbStart = 0.5;
    K_fbEnd   = 0.0;
    K_fbStep  = -0.0009765625;
  }
  
  else {
    cerr << "\nChoose scan type from {scan, upsweep, downsweep}.\n" << endl;
    exit(0);
  }
  
  cout << "phi = " << p.phi << endl;
  
  // coupling strength loop, <= for scan and upsweep, >= for downsweep
  for (K_fb = K_fbStart; K_fb <= K_fbEnd; K_fb += K_fbStep) {
    
    p.K_fb = K_fb;
    
    TS.t = 0.0;
    
    // set initial state of system for new detuning parameter
    if (p.scanType == "scan") TS.N.setInitialConditions();
    
    // time series loop
    for (long long int j = 1; j < p.nSteps; j++) {
      
      // perform time series step
      TS.timeSeriesStep(TS.N, j);
      
      // store time and intensity values
      if (j >= p.transSteps) {
        TS.timeValues(j - p.transSteps)       = TS.t;
        TS.intensityValues0(j - p.transSteps) = norm(TS.N.nodes[0].E);
      }
      
    }
    
    A.fixedPointCheck(TS.intensityValues0);
    A.analyze(TS.timeValues, TS.intensityValues0);
    A.writeFeedbackBifurcationData(file);
    A.resetData();
    
  }
  
  file.close();
  
  // end of time measure
  auto t_end = chrono::steady_clock::now();
  
  cout << "\nFeedback bifurcation line scan was computed with the following parameters:\n" << endl;
  
  TS.showParameters();
  
  cout << "\ncomputation time: " 
       << chrono::duration_cast<chrono::seconds>(t_end - t_start).count() / 60.0
       << " min\n" << endl;
  
}



// generates data for Delta_omega(beta) plot
void rotatingFrameChange ()
{
  
  updateNumericParameters();
  TimeSeries TS;
  
  double betaExp;
  double betaExpStart = -7.0;
  double betaExpEnd   = 0.0;
  double betaExpStep  = 0.0625;
  
  ofstream file("rotating_frame/free_running_laser.txt");
  file.precision(12);
  
  for (betaExp = betaExpStart; betaExp <= betaExpEnd; 
       betaExp += betaExpStep) {
    
    np.beta = pow(10.0, betaExp);
    TS.setRotatingFrame();
    
    file << np.beta << '\t' << np.omega0 << '\t' << np.E0 << endl;
    
  }
  
  file.close();
  
}



// generates a data point for DnuInj_crit(beta) plot
// has to be called with command line options -K_inj and -beta
// use argfileType "SNIPER_change" 
void SNIPER_change (double bifurcationTol, double intensityTol, 
                    bool down)
{
  
  updateNumericParameters();
  TimeSeries TS;
  TISI_Finder T;
  
  double DnuInj_crit;
  string path;
  
  //~ if (p.tau == 0.0) {
    //~ path = "coherence_resonance/DnuInj_crit(beta)_injection.txt";
  //~ }
  
  //~ else if (p.tau > 0.0) {
    //~ path = "coherence_resonance/DnuInj_crit(beta)_injection_feedback.txt";
  //~ }
  
  ofstream file("new.txt");//, ios_base::app);
  file.precision(15);
  
  double betaExpStart     = -7.0;
  double betaExpEnd       = 0.0;
  double betaExpIncrement = 0.0625;
  
  //~ file.open("DnuInj_crit(beta)_injection.txt");
    
  for (double betaExp = betaExpStart; betaExp <= betaExpEnd; 
       betaExp += betaExpIncrement) {
    
    np.beta = pow(10.0, betaExp);
    
    DnuInj_crit = T.findSNIPER(bifurcationTol, intensityTol, down);
    
    cout << "FUCK" << endl;
    cout << "FUCK" << endl;
    cout << "FUCK" << endl;
    cout << "FUCK" << endl;
    cout << "FUCK" << endl;
    cout << "FUCK" << endl;
    cout << "FUCK" << endl;
    cout << "FUCK" << endl;
    cout << "FUCK" << endl;
    cout << "FUCK" << endl;
    
    file << np.beta << '\t' << DnuInj_crit << endl;
    
  }
  
  
  file.close();
  
}



// generates a data point for K_fb_crit(beta) plot
// has to be called with command line options -phi and -beta
// use argfileType "homoclinic_change" 
void homoclinic_change (double bifurcationTol, double intensityTol, 
                        bool left)
{
  
  updateNumericParameters();
  TimeSeries TS;
  TISI_Finder T;
  
  double K_fb_crit;
  
  ofstream file("coherence_resonance/K_fb_crit(beta).txt", 
                ios_base::app);
  file.precision(12);
  
  K_fb_crit = T.findHomoclinic(bifurcationTol, intensityTol, left);
  
  file << np.beta << '\t' << K_fb_crit << endl;
  file.close();
  
}



// generates sampling data for stochastic integration for statistics check
void noiseCheck (int numberOfSamples)
{
  
  complex<double> W;
  
  updateNumericParameters();
  TimeSeries TS;
  
  ofstream file("noise_check/noise_check.txt");
  
  for (int i = 0; i < numberOfSamples; i++) {
    
    W = TS.S.dW(TS.N)(0);
    
    file << W.real() << '\t' << W.imag() << endl;
    cout << i+1 << " samples" << endl;
    
  }
  
  file.close();
  
}



// generates data for one point of injection coherence resonance curve
// has to be called with command line options -K_inj and -beta
// use argfileType "coherenceResonanceInjection"
void TisiVarianceInjection (int numberOfIntervals)
{
  
  // start of time measure
  auto t_start = chrono::steady_clock::now();
  
  updateNumericParameters();
  TimeSeries TS(p.nSteps);
  TISI_Finder T(numberOfIntervals);
  
  int cycleCounter = 0;
  
  int spikeCounter0 = 0;
  int spikeCounter1 = 0;
  int spikeCounter2 = 0;
  int spikeCounter3 = 0;
  
  // define path and file name for coherence resonance file
  string path;
  
  if (p.networkType == "singleLaser" 
   || p.networkType == "singleLaserDelay") {
    
    if (p.tau == 0.0) {
      path = "coherence_resonance/single_laser/injection/coherence_resonance_beta=";
    }
    
    else if (p.tau > 0.0) {
      path = "coherence_resonance/single_laser/injection_feedback/coherence_resonance_beta=";
    }
    
  }
  
  else if (p.networkType == "fourNodeAllToAll") {
    path = "coherence_resonance/network/coherence_resonance_beta=";
  }
  
  path += to_string(np.beta) + ".txt";
  
  ofstream file(path);
  file.precision(12);
  
  cout << "beta = " << np.beta << endl;
  
  // compute critical SNIPER bifurcation parameter
  np.DnuInj = T.findSNIPER();
  
  // set initial state of system
  TS.N.setInitialConditions();
  
  // spike counting loop
  while (true) {
    
    if ((p.networkType == "singleLaser" 
      || p.networkType == "singleLaserDelay") 
      && spikeCounter0 > numberOfIntervals) break;
    
    else if (p.networkType == "fourNodeAllToAll"
          && spikeCounter0 > numberOfIntervals
          && spikeCounter1 > numberOfIntervals
          && spikeCounter2 > numberOfIntervals
          && spikeCounter3 > numberOfIntervals) break;
    
    // stochastic time series loop
    for (long long int i = 0; i < p.nSteps; i++) {
      
      // perform time series step
      TS.timeSeriesStep(TS.N, i + 1 + cycleCounter * p.nSteps);
      
      TS.timeValues(i)       = TS.t;
      TS.intensityValues0(i) = norm(TS.N.nodes[0].E);
      
      if (p.networkType == "fourNodeAllToAll") {
        TS.intensityValues1(i) = norm(TS.N.nodes[1].E);
        TS.intensityValues2(i) = norm(TS.N.nodes[2].E);
        TS.intensityValues3(i) = norm(TS.N.nodes[3].E);
      }
      
    }
    
    if ((p.networkType == "singleLaser" 
     || p.networkType == "singleLaserDelay") && cycleCounter > 0) {
      T.storeSpikeTimes(spikeCounter0, TS.timeValues, TS.intensityValues0);
    }
    
    else if (p.networkType == "fourNodeAllToAll" && cycleCounter > 0) {
      T.storeSpikeTimes(spikeCounter0, spikeCounter1, 
                        spikeCounter2, spikeCounter3, TS.timeValues, 
                        TS.intensityValues0, TS.intensityValues1, 
                        TS.intensityValues2, TS.intensityValues3);
    }
    
    TS.timeValues.zeros();
    TS.intensityValues0.zeros();
    
    if (p.networkType == "fourNodeAllToAll") {
      TS.intensityValues1.zeros();
      TS.intensityValues2.zeros();
      TS.intensityValues3.zeros();
    }
    
    cycleCounter++;
    
    if (p.networkType == "singleLaser" 
     || p.networkType == "singleLaserDelay") {
      cout << p.Ttrans + p.Teval << " ns integrated" << '\t' 
           << cycleCounter       << " times"         << '\t' 
           << spikeCounter0      << " spikes"        << endl;
    }
    
    else if (p.networkType == "fourNodeAllToAll") {
      cout << p.Ttrans + p.Teval << " ns integrated"  << '\t' 
           << cycleCounter       << " times"          << '\t' 
           << spikeCounter0      << " laser 0 spikes" << '\t' 
           << spikeCounter1      << " laser 1 spikes" << '\t' 
           << spikeCounter2      << " laser 2 spikes" << '\t' 
           << spikeCounter3      << " laser 3 spikes" << endl;
    }
    
  }
  
  T.getTISI();
  T.showSpikeTimesAndTISI();
  T.writeCumulantsToFile(file, TS.t);
  
  file.close();
  
  // end of time measure
  auto t_end = chrono::steady_clock::now();
  
  cout << "\nTisi variance for injection was computed with the following parameters:\n" << endl;
  
  TS.showParameters();
  
  cout << "\ncomputation time: " 
       << chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count() / 1e3
       << " s\n" << endl;
  
}



// generates data for one point of feedback coherence resonance curve
// has to be called with command line options -phi and -beta
// use argfileType "coherenceResonanceFeedback"
void TisiVarianceFeedback (int numberOfIntervals) {}
//~ {
  
  //~ // start of time measure
  //~ auto t_start = chrono::steady_clock::now();
  
  //~ updateNumericParameters();
  //~ TimeSeries TS(p.nSteps);
  //~ TISI_Finder T(numberOfIntervals);
  
  //~ int spikeCounter = 0;
  //~ int cycleCounter = 0;
  
  //~ // define path and file name for coherence resonance file
  //~ string path = "coherence_resonance/coherence_resonance_feedback_beta=" 
              //~ + to_string(np.beta) + ".txt";
  //~ ofstream file(path);
  //~ file.precision(12);
  
  //~ cout << "beta = " << np.beta << endl;
  
  //~ // compute critical homoclinic bifurcation parameter
  //~ p.K_fb = T.findHomoclinic();
  
  //~ p.timeSeriesMode    = "stochastic";
  //~ p.integrationMethod = "explicitEulerMaruyamaDelay";
  //~ p.initialCondition  = "small";
  //~ p.E_sp_det          = false;
  
  //~ // set initial state of system
  //~ TS.N.setInitialConditions();
  
  //~ // spike counting loop
  //~ while (spikeCounter < numberOfIntervals && TS.t < 1e7) {
    
    //~ // loop for stochastic time series
    //~ for (long long int i = 0; i < p.nSteps; i++) {
      
      //~ // perform time series step
      //~ TS.timeSeriesStep(TS.N, i + 1 + cycleCounter * p.nSteps);
      
      //~ TS.timeValues(i)       = TS.t;
      //~ TS.intensityValues0(i) = norm(TS.N.nodes[0].E);
      
    //~ }
    
    //~ if (cycleCounter > 0) {
      //~ T.storeSpikeTimes(spikeCounter, TS.timeValues, TS.intensityValues0);
    //~ }
    
    //~ TS.timeValues.zeros();
    //~ TS.intensityValues0.zeros();
    
    //~ cycleCounter++;
    
    //~ cout << p.Ttrans + p.Teval << " ns integrated"   << '\t' 
         //~ << cycleCounter       << " times"           << '\t' 
         //~ << spikeCounter       << " spikes detected" << endl;
    
  //~ }
  
  //~ T.getTISI();
  //~ T.showSpikeTimesAndTISI();
  //~ T.writeCumulantsToFile(file, TS.t);
  
  //~ file.close();
  
  //~ // end of time measure
  //~ auto t_end = chrono::steady_clock::now();
  
  //~ cout << "\nTisi variance for feedback was computed with the following parameters:\n" << endl;
  
  //~ TS.showParameters();
  
  //~ cout << "\ncomputation time: " 
       //~ << chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count() / 1e3
       //~ << " s\n" << endl;
  
//~ }

