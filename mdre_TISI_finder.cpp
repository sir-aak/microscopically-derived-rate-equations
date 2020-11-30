#include "mdre_parameters.hpp"
#include "mdre_network.hpp"
#include "mdre_solver.hpp"
#include "mdre_auxiliary.hpp"
#include "mdre_time_series.hpp"
#include "mdre_TISI_finder.hpp"


extern NodeParameters np;
extern Parameters p;


// TISI_finder constructor
TISI_Finder::TISI_Finder ()
{
  
  numberOfIntervals  = 1e3;
  
  spikeThreshold_up  = 1.60;
  spikeThreshold_low = 0.70;
  
  tlimit0 = tlimit1 = tlimit2 = tlimit3 = 0.0;
  tdead   = 0.25;
  
  // simulates dead-time after spike detection
  wait0 = wait1 = wait2 = wait3 = false;
  
  spikeTimes0.zeros(numberOfIntervals + 1);
  Tisi0.zeros(numberOfIntervals);
  
  if (p.networkType == "fourNodeAllToAll") {
    
    spikeTimes1.zeros(numberOfIntervals + 1);
    Tisi1.zeros(numberOfIntervals);
    
    spikeTimes2.zeros(numberOfIntervals + 1);
    Tisi2.zeros(numberOfIntervals);
    
    spikeTimes3.zeros(numberOfIntervals + 1);
    Tisi3.zeros(numberOfIntervals);
    
  }
  
}



// overloaded input-version
TISI_Finder::TISI_Finder (int _numberOfIntervals)
{
  
  numberOfIntervals  = _numberOfIntervals;
  
  spikeThreshold_up  = 1.60;
  spikeThreshold_low = 0.70;
  
  tlimit0 = tlimit1 = tlimit2 = tlimit3 = 0.0;
  tdead  = 0.25;
  
  // simulates dead-time after spike detection
  wait0 = wait1 = wait2 = wait3 = false;
  
  spikeTimes0.zeros(numberOfIntervals + 1);
  Tisi0.zeros(numberOfIntervals);
  
  if (p.networkType == "fourNodeAllToAll") {
    
    spikeTimes1.zeros(numberOfIntervals + 1);
    Tisi1.zeros(numberOfIntervals);
    
    spikeTimes2.zeros(numberOfIntervals + 1);
    Tisi2.zeros(numberOfIntervals);
    
    spikeTimes3.zeros(numberOfIntervals + 1);
    Tisi3.zeros(numberOfIntervals);
    
  }
  
}



// returns down-rounded value of a to the decimal's decimal place
double TISI_Finder::roundDown (double a, int decimal)
{
  
  double b = pow(10, decimal);
  
  if (a > 0.0) return (floor(a * b) / b);
  else         return (ceil(a * b) / b);
  
}



// returns up-rounded value of a to the decimal's decimal place
double TISI_Finder::roundUp (double a, int decimal)
{
  
  double b = pow(10, decimal);
  
  if (a > 0.0) return (ceil(a * b) / b);
  else         return (floor(a * b) / b);
  
}



// returns DnuInj_crit where SNIPER bifurcation occurs at given K_inj and beta
// has to be called with the command line options -K_inj and -beta
double TISI_Finder::findSNIPER (double bifurcationTol, 
                                double intensityTol, bool down)
{
  
  // start of time measure
  auto t_start = chrono::steady_clock::now();
  
  if (np.K_inj == 0.0) {
    cerr << "\nThere is no SNIPER bifurcation for K_inj = 0.\n" << endl;
    exit(0);
  }
  
  // store original parameter values and modes
  string timeSeriesMode    = p.timeSeriesMode;
  string integrationMethod = p.integrationMethod;
  bool   E_sp_det          = p.E_sp_det;
  double dt                = p.dt;
  
  
  // set parameters for computation of critical bifurcation parameter
  
  if (p.tau == 0.0) {
    p.timeSeriesMode    = "deterministic";
    p.integrationMethod = "explicitRungeKutta4";
    p.E_sp_det          = true;
    p.dt                = 1e-3;
  }
  
  else if (p.tau > 0.0) {
    p.timeSeriesMode    = "deterministic";
    p.integrationMethod = "explicitRungeKutta4Delay";
    p.E_sp_det          = true;
    p.dt                = 1e-3;
  }
  
  updateNumericParameters();
  TimeSeries TS;
  
  
  // precision for output
  int precision = int(-log10(bifurcationTol)) + 2;
  
  // for additional check of dynamics type
  bool   earlySpike          = false;
  double earlySpikeThreshold = 1.5;
  int    checkSteps          = int(5.0 / p.dt);
  
  double DnuInj_prev = -1.0;	// start point
  double DnuInjStep;
  
  if (down == true) DnuInjStep = -1.0;
  else              DnuInjStep =  1.0;
  
  np.DnuInj = DnuInj_prev + DnuInjStep;
  
  int iterationCounter = 0;
  
  // loop for finding critical SNIPER bifurcation parameter DnuInj_crit
  while (abs(np.DnuInj - DnuInj_prev) > bifurcationTol / 10.0) {
    
    iterationCounter++;
    DnuInj_prev = np.DnuInj;
    
    cout << "iteration " << iterationCounter        << '\t' 
         << "DnuInj = "  << setprecision(precision) << np.DnuInj << '\t';
    
    // set initial state of system
    TS.t = 0.0;
    TS.N.setInitialConditions();
    
    // time series loop
    for (long long int i = 1; i < p.nSteps; i++) {
      
      // perform time series step
      TS.timeSeriesStep(TS.N, i);
      
      // pertubation to check if solution stays synchronous
      if (TS.t >= 50.0 && TS.t < 55.0 
          && p.networkType == "fourNodeAllToAll" 
          && p.initialCondition == "synchron") {
        TS.networkPerturbation();
      }
      
      // check if early spike occurs when t is in [5 ns, Ttrans]
      if (i >= checkSteps && i <= p.transSteps 
          && norm(TS.N.nodes[0].E) > earlySpikeThreshold * norm(np.E0)) {
        earlySpike = true;
        break;
      }
      
      // store time and intensity values after transient time
      if (i >= p.transSteps) {
        TS.timeValues(i - p.transSteps)       = TS.t;
        TS.intensityValues0(i - p.transSteps) = norm(TS.N.nodes[0].E);
      }
      
    }
    
    // check for fixed point solution
    if (abs(TS.intensityValues0.max() - TS.intensityValues0.min()) 
        < intensityTol && earlySpike == false) {
      
      np.DnuInj += DnuInjStep;
      
      cout << "fixed point" << '\t' << "max - min = " 
           << abs(TS.intensityValues0.max() - TS.intensityValues0.min()) 
           << '\t' << "early spike: " << earlySpike << endl;
      
    }
    
    else {
      
      np.DnuInj  = np.DnuInj - DnuInjStep;
      DnuInjStep = DnuInjStep / 4.0;
      
      cout << "spiking sol" << '\t' << "max - min = " 
           << abs(TS.intensityValues0.max() - TS.intensityValues0.min()) << '\t' 
           << "early spike: " << earlySpike << endl;
      
    }
    
    earlySpike = false;
    
  }
  
  // end of time measure
  auto t_end = chrono::steady_clock::now();
  
  cout << "\nSNIPER bifurcation parameter DnuInj_SNIPER = " 
       << setprecision(precision) 
       << roundDown(DnuInj_prev, precision - 2) << endl;
  cout << "was computed for a tolerance of "    << bifurcationTol << endl;
  cout << "with the following parameters:\n"      << endl;
  
  TS.showParameters();
  
  cout << "\ncomputation time: " 
       << chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count() / 1e3 
       << " s\n" << endl;
  
  // reset parameters to original values
  //~ p.networkType       = networkType;
  p.timeSeriesMode    = timeSeriesMode;
  p.integrationMethod = integrationMethod;
  //~ p.initialCondition  = initialCondition;
  p.E_sp_det          = E_sp_det;
  p.dt                = dt;
  
  return (roundDown(DnuInj_prev, precision - 2));
  
}



// returns K_fb_crit where homoclinic bifurcation occurs at given phi
// and beta
// has to be called with the command line options -phi and -beta
double TISI_Finder::findHomoclinic (double bifurcationTol, 
                                    double intensityTol, bool left)
{
  
  // start of time measure
  auto t_start = chrono::steady_clock::now();
  
  // store original parameter values and modes
  string networkType       = p.networkType;
  string timeSeriesMode    = p.timeSeriesMode;
  string integrationMethod = p.integrationMethod;
  string initialCondition  = p.initialCondition;
  bool   E_sp_det          = p.E_sp_det;
  
  
  // set parameters for computation of critical bifurcation parameter
  
  if (p.tau == 0.0) {
    p.networkType       = "singleLaser";
    p.timeSeriesMode    = "deterministic";
    p.integrationMethod = "explicitRungeKutta4";
    p.initialCondition  = "small";
    p.E_sp_det          = true;
  }
  
  else if (p.tau > 0.0) {
    p.networkType       = "singleLaserDelay";
    p.timeSeriesMode    = "deterministic";
    p.integrationMethod = "explicitRungeKutta4Delay";
    p.initialCondition  = "small";
    p.E_sp_det          = true;
  }
  
  updateNumericParameters();
  TimeSeries TS;
  
  
  // precision for output
  int precision = int(-log10(bifurcationTol)) + 2;
  
  double K_fb_prev = 0.15;
  double K_fb_Step;
  
  if (left == true) K_fb_Step = -0.01;
  else              K_fb_Step =  0.01;
  
  p.K_fb = K_fb_prev + K_fb_Step;
  
  int iterationCounter = 0;
  
  // loop for finding critical SNIPER bifurcation parameter DnuInj_crit
  while (abs(p.K_fb - K_fb_prev) > bifurcationTol / 10.0) {
    
    iterationCounter++;
    K_fb_prev = p.K_fb;
    
    cout << "iteration " << iterationCounter        << '\t' 
         << "K_fb = "    << setprecision(precision) << p.K_fb << '\t';
    
    // set initial state of system
    TS.t = 0.0;
    TS.N.setInitialConditions();
    
    // time series loop
    for (long long int i = 1; i < p.nSteps; i++) {
      
      // perform time series step
      TS.timeSeriesStep(TS.N, i);
      
      // store time and intensity values after transient time
      if (i >= p.transSteps) {
        TS.timeValues(i - p.transSteps)       = TS.t;
        TS.intensityValues0(i - p.transSteps) = norm(TS.N.nodes[0].E);
      }
      
    }
    
    // check for fixed point solution
    if (abs(TS.intensityValues0.max() - TS.intensityValues0.min()) 
        < intensityTol) { // && earlySpike == false) {
      
      p.K_fb += K_fb_Step;
      
      cout << "fixed point" << '\t' << "max - min = " 
           << abs(TS.intensityValues0.max() - TS.intensityValues0.min()) << endl;
      
    }
    
    else {
      
      p.K_fb    = p.K_fb - K_fb_Step;
      K_fb_Step = K_fb_Step / 4.0;
      
      cout << "spiking sol" << '\t' << "max - min = " 
           << abs(TS.intensityValues0.max() - TS.intensityValues0.min()) << endl;
      
    }
    
  }
  
  // end of time measure
  auto t_end = chrono::steady_clock::now();
  
  cout << "\nHomoclinic bifurcation parameter K_fb_homoclinic = " 
       << setprecision(precision) 
       << roundDown(K_fb_prev, precision - 2) << endl;
  cout << "was computed for a tolerance of "  << bifurcationTol << endl;
  cout << "with the following parameters:\n"  << endl;
  
  TS.showParameters();
  
  cout << "\ncomputation time: " 
       << chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count() / 1e3 
       << " s\n" << endl;
  
  // reset parameters to original values
  p.networkType       = networkType;
  p.timeSeriesMode    = timeSeriesMode;
  p.integrationMethod = integrationMethod;
  p.initialCondition  = initialCondition;
  p.E_sp_det          = E_sp_det;
  
  return (roundUp(K_fb_prev, precision - 2));
  
}



// extracts all spike times from single laser time series
void TISI_Finder::storeSpikeTimes (int& spikeCounter0, 
                                   arma::colvec& timeValues, 
                                   arma::colvec& intensityValues0)
{
  
  // in time series, check if ..
  for (int i = 0; i < p.nSteps; i++) {
    
    // .. enough spikes are detected
    if (spikeCounter0 == numberOfIntervals + 1) break;
    
    // .. there is a spike in laser
    // and if yes, store corresponding time value
    if (intensityValues0(i) < spikeThreshold_low 
     && timeValues(i) > tlimit0 
     && spikeCounter0 < numberOfIntervals + 1) {
      tlimit0 = timeValues(i) + tdead;
      spikeTimes0(spikeCounter0) = timeValues(i);
      spikeCounter0++;
    }
    
  }
  
}



// extracts all spike times from laser network time series
void TISI_Finder::storeSpikeTimes (int& spikeCounter0, int& spikeCounter1, 
                                   int& spikeCounter2, int& spikeCounter3, 
                                   arma::colvec& timeValues, 
                                   arma::colvec& intensityValues0, 
                                   arma::colvec& intensityValues1, 
                                   arma::colvec& intensityValues2, 
                                   arma::colvec& intensityValues3)
{
  
  // in time series, check if ..
  for (int i = 0; i < p.nSteps; i++) {
    
    // .. enough spikes are detected
    if (spikeCounter0 == numberOfIntervals + 1 
     && spikeCounter1 == numberOfIntervals + 1 
     && spikeCounter2 == numberOfIntervals + 1 
     && spikeCounter3 == numberOfIntervals + 1) break;
    
    
    // laser 0
    
    if (spikeCounter0 < numberOfIntervals + 1) {
      
      // .. there is a spike in laser 0
      // and if yes, store corresponding time value
      if (intensityValues0(i) < spikeThreshold_low 
       && timeValues(i) > tlimit0 
       && spikeCounter0 < numberOfIntervals + 1) {
        tlimit0 = timeValues(i) + tdead;
        spikeTimes0(spikeCounter0) = timeValues(i);
        spikeCounter0++;
      }
      
    }
    
    // *****************************************************************
    
    // laser 1
    
    if (spikeCounter1 < numberOfIntervals + 1) {
      
      // .. there is a spike in laser 1
      // and if yes, store corresponding time value
      if (intensityValues1(i) < spikeThreshold_low 
       && timeValues(i) > tlimit1 
       && spikeCounter1 < numberOfIntervals + 1) {
        tlimit1 = timeValues(i) + tdead;
        spikeTimes1(spikeCounter1) = timeValues(i);
        spikeCounter1++;
      }
      
    }
    
    // *****************************************************************
    
    // laser 2
    
    if (spikeCounter2 < numberOfIntervals + 1) {
      
      // .. there is a spike in laser 2
      // and if yes, store corresponding time value
      if (intensityValues2(i) < spikeThreshold_low 
       && timeValues(i) > tlimit2 
       && spikeCounter2 < numberOfIntervals + 1) {
        tlimit2 = timeValues(i) + tdead;
        spikeTimes2(spikeCounter2) = timeValues(i);
        spikeCounter2++;
      }
      
    }
    
    // *****************************************************************
    
    // laser 3
    
    if (spikeCounter3 < numberOfIntervals + 1) {
      
      // .. there is a spike in laser 3
      // and if yes, store corresponding time value
      if (intensityValues3(i) < spikeThreshold_low 
       && timeValues(i) > tlimit3 
       && spikeCounter3 < numberOfIntervals + 1) {
        tlimit3 = timeValues(i) + tdead;
        spikeTimes3(spikeCounter3) = timeValues(i);
        spikeCounter3++;
      }
      
    }
    
  }
  
}



// calculates inter-spike intervals
void TISI_Finder::getTISI ()
{
  
  Tisi0 = diff(spikeTimes0);
  
  if (p.networkType == "fourNodeAllToAll") {
    Tisi1 = diff(spikeTimes1);
    Tisi2 = diff(spikeTimes2);
    Tisi3 = diff(spikeTimes3);
  }
  
}



// shows all spike times and inter-spike intervals
void TISI_Finder::showSpikeTimesAndTISI ()
{
  
  if (p.networkType == "singleLaser" 
   || p.networkType == "singleLaserDelay") {
    
    cout << "\nspike times:\n" << endl;
    
    for (unsigned int i = 0; i < spikeTimes0.n_elem; i++) {
      cout << spikeTimes0(i) << endl;
    }
    
    cout << "\nTisi:\n" << endl;
    
    for (unsigned int i = 0; i < Tisi0.n_elem; i++) {
      cout << Tisi0(i) << endl;
    }
    
  }
  
  else if (p.networkType == "fourNodeAllToAll") {
    
    cout << "\nspike times:\n" << endl;
    
    for (unsigned int i = 0; i < spikeTimes0.n_elem; i++) {
      cout << spikeTimes0(i) << '\t' << spikeTimes1(i) << '\t' 
           << spikeTimes2(i) << '\t' << spikeTimes3(i) << endl;
    }
    
    cout << "\nTisi:\n" << endl;
    
    for (unsigned int i = 0; i < Tisi0.n_elem; i++) {
      cout << Tisi0(i) << '\t' << Tisi1(i) << '\t' 
           << Tisi2(i) << '\t' << Tisi3(i) << endl;
    }
    
  }
  
}



// writes mean, variance and all inter-spike intervals in file
void TISI_Finder::writeCumulantsToFile (ofstream& file, double t)
{
  
  if (p.networkType == "singleLaser" 
   || p.networkType == "singleLaserDelay") {
    
    file << np.beta     << '\t' 
         << mean(Tisi0) << '\t' << var(Tisi0, 1) << '\t' 
         << t << '\n'   << endl;
    
    for (unsigned int i = 0; i < Tisi0.n_elem; i++) {
      file << Tisi0(i) << endl;
    }
    
  }
  
  else if (p.networkType == "fourNodeAllToAll") {
    
    file << np.beta     << '\t' 
         << mean(Tisi0) << '\t' << var(Tisi0, 1) << '\t' 
         << mean(Tisi1) << '\t' << var(Tisi1, 1) << '\t' 
         << mean(Tisi2) << '\t' << var(Tisi2, 1) << '\t' 
         << mean(Tisi3) << '\t' << var(Tisi3, 1) << '\t' 
         << t << '\n'   << endl;
    
    for (unsigned int i = 0; i < Tisi0.n_elem; i++) {
      file << Tisi0(i) << '\t' << Tisi1(i) << '\t' 
           << Tisi2(i) << '\t' << Tisi3(i) << endl;
    }
    
  }
  
}

