#include "mdre_parameters.hpp"
#include "mdre_network.hpp"
#include "mdre_analyzer.hpp"


extern NodeParameters np;
extern Parameters p;


// constructor for Analyzer class
// purpose of this class is: 
// check for fixed point state and get cwIntensity
// finding and counting of extrema and determine period of time series
Analyzer::Analyzer ()
{
  resetData();
}



// sets all class parameters back to zero, empties extrema lists
void Analyzer::resetData ()
{
  
  fixedPoint0 = false;
  fixedPoint1 = false;
  fixedPoint2 = false;
  fixedPoint3 = false;
  
  cwIntensity0 = NAN;
  cwIntensity1 = NAN;
  cwIntensity2 = NAN;
  cwIntensity3 = NAN;
  
  maximaCounter0 = 0;       minimaCounter0 = 0;
  maximaCounter1 = 0;       minimaCounter1 = 0;
  maximaCounter2 = 0;       minimaCounter2 = 0;
  maximaCounter3 = 0;       minimaCounter3 = 0;
  
  uniqueMaximaCounter0 = 0; uniqueMinimaCounter0 = 0;
  uniqueMaximaCounter1 = 0; uniqueMinimaCounter1 = 0;
  uniqueMaximaCounter2 = 0; uniqueMinimaCounter2 = 0;
  uniqueMaximaCounter3 = 0; uniqueMinimaCounter3 = 0;
  
  maxima0.clear();          maximaTime0.clear();
  maxima1.clear();          maximaTime1.clear();
  maxima2.clear();          maximaTime2.clear();
  maxima3.clear();          maximaTime3.clear();
  
  minima0.clear();          minimaTime0.clear();
  minima1.clear();          minimaTime1.clear();
  minima2.clear();          minimaTime2.clear();
  minima3.clear();          minimaTime3.clear();
  
  uniqueMaxima0.clear();    uniqueMinima0.clear();
  uniqueMaxima1.clear();    uniqueMinima1.clear();
  uniqueMaxima2.clear();    uniqueMinima2.clear();
  uniqueMaxima3.clear();    uniqueMinima3.clear();
  
  period0 = 0.0;
  period1 = 0.0;
  period2 = 0.0;
  period3 = 0.0;
  
  deviation01 = 0.0;
  deviation02 = 0.0;
  deviation03 = 0.0;
  deviation12 = 0.0;
  deviation13 = 0.0;
  deviation23 = 0.0;
  
  LyapunovExponent = NAN;
  
}



// checks single laser time series for fixed point solution
void Analyzer::fixedPointCheck (arma::colvec& intensityValues0, 
                                double tol)
{
  
  if (intensityValues0.max() - intensityValues0.min() < tol) {
    fixedPoint0          = true;
    uniqueMaximaCounter0 = 0;
    uniqueMinimaCounter0 = 0;
    period0              = 0.0;
  }
  
}



// checks network time series for fixed point solutions
void Analyzer::fixedPointCheck (arma::colvec& intensityValues0, 
                                arma::colvec& intensityValues1, 
                                arma::colvec& intensityValues2, 
                                arma::colvec& intensityValues3, 
                                double tol)
{
  
  if (intensityValues0.max() - intensityValues0.min() < tol) {
    fixedPoint0          = true;
    uniqueMaximaCounter0 = 0;
    uniqueMinimaCounter0 = 0;
    period0              = 0.0;
  }
  
  if (intensityValues1.max() - intensityValues1.min() < tol) {
    fixedPoint1          = true;
    uniqueMaximaCounter1 = 0;
    uniqueMinimaCounter1 = 0;
    period1              = 0.0;
  }
  
  if (intensityValues2.max() - intensityValues2.min() < tol) {
    fixedPoint2          = true;
    uniqueMaximaCounter2 = 0;
    uniqueMinimaCounter2 = 0;
    period2              = 0.0;
  }
  
  if (intensityValues3.max() - intensityValues3.min() < tol) {
    fixedPoint3          = true;
    uniqueMaximaCounter3 = 0;
    uniqueMinimaCounter3 = 0;
    period3              = 0.0;
  }
  
}



// sets laser intensity of fixed point solution for single laser
void Analyzer::setcwIntensity (arma::colvec& intensityValues0)
{
  if (fixedPoint0 == true) cwIntensity0 = intensityValues0.back();
  else                     cwIntensity0 = NAN;
}



// sets laser intensity of fixed point solution for laser network
void Analyzer::setcwIntensity (arma::colvec& intensityValues0, 
                               arma::colvec& intensityValues1, 
                               arma::colvec& intensityValues2, 
                               arma::colvec& intensityValues3)
{
  
  if (fixedPoint0 == true) cwIntensity0 = intensityValues0.back();
  else                     cwIntensity0 = NAN;
  
  if (fixedPoint1 == true) cwIntensity1 = intensityValues1.back();
  else                     cwIntensity1 = NAN;
  
  if (fixedPoint2 == true) cwIntensity2 = intensityValues2.back();
  else                     cwIntensity2 = NAN;
  
  if (fixedPoint3 == true) cwIntensity3 = intensityValues3.back();
  else                     cwIntensity3 = NAN;
  
}



// returns "true" if a maximum is found in time series
bool Analyzer::isMaximum (double I_prev, double I_curr, double I_next, 
                          double tol)
{
  return (I_curr > I_prev && I_curr > I_next 
       && abs(I_curr - I_prev) > tol && abs(I_curr - I_next) > tol);
}



// returns "true" if a minimum is found in time series
bool Analyzer::isMinimum (double I_prev, double I_curr, double I_next, 
                          double tol)
{
  return (I_curr < I_prev && I_curr < I_next 
       && abs(I_curr - I_prev) > tol && abs(I_curr - I_next) > tol);
}



// returns "true" if intensity value I already exists in extremaValues-list
bool Analyzer::extremumExists (double I, vector<double> uniqueExtrema, 
                               double tol)
{
  
  if (uniqueExtrema.empty()) return (false);
  
  else {
    
    for (unsigned int k = 0; k < uniqueExtrema.size(); k++) {
      
      if (abs(I - uniqueExtrema[k]) < tol) {
        return (true);
      }
      
    }
    
  }
  
  return (false);
  
}



// returns period from given maxima and maximaTime lists
double Analyzer::getPeriod (vector<double> maxima, vector<double> maximaTime, 
                            double tol)
{
  
  double period = 0.0;
  
  for (unsigned int i = 1; i < maxima.size(); i++) {
    
    if ((maxima[i] <= maxima[0] + tol) 
     && (maxima[i] >= maxima[0] - tol)) {
      period = maximaTime[i] - maximaTime[0];
      break;
    }
    
  }
  
  return (period);
  
}



// counts and stores extrema for single laser
void Analyzer::countAndStoreExtrema (int k, arma::colvec& timeValues, 
                                     arma::colvec& intensityValues0)
{
  
  if (fixedPoint0 == false) {
    
    // .. maximum is found
    if (isMaximum(intensityValues0(k-1), intensityValues0(k), 
                  intensityValues0(k+1))) {
      maxima0.push_back(intensityValues0(k));
      maximaTime0.push_back(timeValues(k));
      maximaCounter0++;
    }
    
    // .. minimum is found
    if (isMinimum(intensityValues0(k-1), intensityValues0(k), 
                  intensityValues0(k+1))) {
      minima0.push_back(intensityValues0(k));
      minimaTime0.push_back(timeValues(k));
      minimaCounter0++;
    }
    
    // .. unique maximum is found
    if (isMaximum(intensityValues0(k-1), intensityValues0(k), 
                  intensityValues0(k+1)) 
     && !extremumExists(intensityValues0(k), uniqueMaxima0)) {
      uniqueMaxima0.push_back(intensityValues0(k));
      uniqueMaximaCounter0++;
    }
    
    // .. unique minimum is found
    if (isMinimum(intensityValues0(k-1), intensityValues0(k), 
                  intensityValues0(k+1)) 
     && !extremumExists(intensityValues0(k), uniqueMinima0)) {
      uniqueMinima0.push_back(intensityValues0(k));
      uniqueMinimaCounter0++;
    }
    
  }
  
}



// counts and stores extrema for laser network
void Analyzer::countAndStoreExtrema (int k, arma::colvec& timeValues, 
     arma::colvec& intensityValues0, arma::colvec& intensityValues1, 
     arma::colvec& intensityValues2, arma::colvec& intensityValues3)
{
  
  // checking laser 0
  if (fixedPoint0 == false) {
    
    // .. maximum is found
    if (isMaximum(intensityValues0(k-1), intensityValues0(k), 
                  intensityValues0(k+1))) {
      maxima0.push_back(intensityValues0(k));
      maximaTime0.push_back(timeValues(k));
      maximaCounter0++;
    }
    
    // .. minimum is found
    if (isMinimum(intensityValues0(k-1), intensityValues0(k), 
                  intensityValues0(k+1))) {
      minima0.push_back(intensityValues0(k));
      minimaTime0.push_back(timeValues(k));
      minimaCounter0++;
    }
    
    // .. unique maximum is found
    if (isMaximum(intensityValues0(k-1), intensityValues0(k), 
                  intensityValues0(k+1)) 
     && !extremumExists(intensityValues0(k), uniqueMaxima0)) {
      uniqueMaxima0.push_back(intensityValues0(k));
      uniqueMaximaCounter0++;
    }
    
    // .. unique minimum is found
    if (isMinimum(intensityValues0(k-1), intensityValues0(k), 
                  intensityValues0(k+1)) 
     && !extremumExists(intensityValues0(k), uniqueMinima0)) {
      uniqueMinima0.push_back(intensityValues0(k));
      uniqueMinimaCounter0++;
    }
    
  }
  
  // *******************************************************************
  
  // checking laser 1
  if (fixedPoint1 == false) {
    
    // .. maximum is found
    if (isMaximum(intensityValues1(k-1), intensityValues1(k), 
                  intensityValues1(k+1))) {
      maxima1.push_back(intensityValues1(k));
      maximaTime1.push_back(timeValues(k));
      maximaCounter1++;
    }
    
    // .. minimum is found
    if (isMinimum(intensityValues1(k-1), intensityValues1(k), 
                  intensityValues1(k+1))) {
      minima1.push_back(intensityValues1(k));
      minimaTime1.push_back(timeValues(k));
      minimaCounter1++;
    }
    
    // .. unique maximum is found
    if (isMaximum(intensityValues1(k-1), intensityValues1(k), 
                  intensityValues1(k+1)) 
     && !extremumExists(intensityValues1(k), uniqueMaxima1)) {
      uniqueMaxima1.push_back(intensityValues1(k));
      uniqueMaximaCounter1++;
    }
    
    // .. unique minimum is found
    if (isMinimum(intensityValues1(k-1), intensityValues1(k), 
                  intensityValues1(k+1)) 
     && !extremumExists(intensityValues1(k), uniqueMinima1)) {
      uniqueMinima1.push_back(intensityValues1(k));
      uniqueMinimaCounter1++;
    }
    
  }
  
  // *******************************************************************
  
  // checking laser 2
  if (fixedPoint2 == false) {
    
    // .. maximum is found
    if (isMaximum(intensityValues2(k-1), intensityValues2(k), 
                  intensityValues2(k+1))) {
      maxima2.push_back(intensityValues2(k));
      maximaTime2.push_back(timeValues(k));
      maximaCounter2++;
    }
    
    // .. minimum is found
    if (isMinimum(intensityValues2(k-1), intensityValues2(k), 
                  intensityValues2(k+1))) {
      minima2.push_back(intensityValues2(k));
      minimaTime2.push_back(timeValues(k));
      minimaCounter2++;
    }
    
    // .. unique maximum is found
    if (isMaximum(intensityValues2(k-1), intensityValues2(k), 
                  intensityValues2(k+1)) 
     && !extremumExists(intensityValues2(k), uniqueMaxima2)) {
      uniqueMaxima2.push_back(intensityValues2(k));
      uniqueMaximaCounter2++;
    }
    
    // .. unique minimum is found
    if (isMinimum(intensityValues2(k-1), intensityValues2(k), 
                  intensityValues2(k+1)) 
     && !extremumExists(intensityValues2(k), uniqueMinima2)) {
      uniqueMinima2.push_back(intensityValues2(k));
      uniqueMinimaCounter2++;
    }
    
  }
  
  // *******************************************************************
  
  // checking laser 3
  if (fixedPoint3 == false) {
    
    // .. maximum is found
    if (isMaximum(intensityValues3(k-1), intensityValues3(k), 
                  intensityValues3(k+1))) {
      maxima3.push_back(intensityValues3(k));
      maximaTime3.push_back(timeValues(k));
      maximaCounter3++;
    }
    
    // .. minimum is found
    if (isMinimum(intensityValues3(k-1), intensityValues3(k), 
                  intensityValues3(k+1))) {
      minima3.push_back(intensityValues3(k));
      minimaTime3.push_back(timeValues(k));
      minimaCounter3++;
    }
    
    // .. unique maximum is found
    if (isMaximum(intensityValues3(k-1), intensityValues3(k), 
                  intensityValues3(k+1)) 
     && !extremumExists(intensityValues3(k), uniqueMaxima3)) {
      uniqueMaxima3.push_back(intensityValues3(k));
      uniqueMaximaCounter3++;
    }
    
    // .. unique minimum is found
    if (isMinimum(intensityValues3(k-1), intensityValues3(k), 
                  intensityValues3(k+1)) 
     && !extremumExists(intensityValues3(k), uniqueMinima3)) {
      uniqueMinima3.push_back(intensityValues3(k));
      uniqueMinimaCounter3++;
    }
    
  }
  
}



// analyzes time series for single laser
void Analyzer::analyze (arma::colvec& timeValues, 
                        arma::colvec& intensityValues0)
{
  
  setcwIntensity(intensityValues0);
  
  for (int k = 1; k < p.evalSteps - 1; k++) {
    countAndStoreExtrema(k, timeValues, intensityValues0);
  }
  
  // calculate period
  period0 = getPeriod(maxima0, maximaTime0);
  
}



// analyzes time series for single laser with Lyapunov exponent
void Analyzer::analyze (arma::colvec& timeValues, 
                        arma::colvec& intensityValues0, 
                        arma::colvec& LyapunovValues)
{
  
  setcwIntensity(intensityValues0);
  
  for (int k = 1; k < p.evalSteps - 1; k++) {
    countAndStoreExtrema(k, timeValues, intensityValues0);
  }
  
  // calculate period
  period0 = getPeriod(maxima0, maximaTime0);
  
  if (fixedPoint0 == true) LyapunovExponent = LyapunovValues.back();
  else                     LyapunovExponent = max(LyapunovValues);
  
}



// analyzes time series for laser network
void Analyzer::analyze (arma::colvec& timeValues, 
     arma::colvec& intensityValues0, arma::colvec& intensityValues1, 
     arma::colvec& intensityValues2, arma::colvec& intensityValues3)
{
  
  setcwIntensity(intensityValues0, intensityValues1, 
                 intensityValues2, intensityValues3);
  
  for (int k = 1; k < p.evalSteps - 1; k++) {
    countAndStoreExtrema(k, timeValues, 
                         intensityValues0, intensityValues1, 
                         intensityValues2, intensityValues3);
  }
  
  
  // calculate period
  
  period0 = getPeriod(maxima0, maximaTime0);
  period1 = getPeriod(maxima1, maximaTime1);
  period2 = getPeriod(maxima2, maximaTime2);
  period3 = getPeriod(maxima3, maximaTime3);
  
  
  // calculate synchronization
  
  deviation01 = norm(intensityValues0 - intensityValues1) 
              / sqrt(intensityValues0.n_elem);
  deviation02 = norm(intensityValues0 - intensityValues2) 
              / sqrt(intensityValues0.n_elem);
  deviation03 = norm(intensityValues0 - intensityValues3) 
              / sqrt(intensityValues0.n_elem);
  deviation12 = norm(intensityValues1 - intensityValues2) 
              / sqrt(intensityValues0.n_elem);
  deviation13 = norm(intensityValues1 - intensityValues3) 
              / sqrt(intensityValues0.n_elem);
  deviation23 = norm(intensityValues2 - intensityValues3) 
              / sqrt(intensityValues0.n_elem);
  
}



void Analyzer::writeInjectionSingleLaserData (ofstream& file)
{
  
  file << np.K_inj                  << '\t' 
       << np.DnuInj                 << '\t' 
       << cwIntensity0              << '\t' 
       << uniqueMaximaCounter0      << '\t' 
       << uniqueMinimaCounter0      << '\t' 
       << period0                   << '\t' 
       << LyapunovExponent          << endl;
  
  cout << np.K_inj                  << '\t' 
       << np.DnuInj                 << '\t' 
       << cwIntensity0              << '\t' 
       << uniqueMaximaCounter0      << '\t' 
       << uniqueMinimaCounter0      << '\t' 
       << period0                   << '\t' 
       << LyapunovExponent          << endl;
  
}



void Analyzer::writeInjectionLaserNetworkData (ofstream& file)
{
  
  file << np.K_inj                  << '\t' // [0]
       << np.DnuInj                 << '\t' // [1]
       << cwIntensity0              << '\t' // [2]
       << uniqueMaximaCounter0      << '\t' // [3]
       << uniqueMinimaCounter0      << '\t' // [4]
       << period0                   << '\t' // [5]
       << cwIntensity1              << '\t' // [6]
       << uniqueMaximaCounter1      << '\t' // [7]
       << uniqueMinimaCounter1      << '\t' // [8]
       << period1                   << '\t' // [9]
       << cwIntensity2              << '\t' // [10]
       << uniqueMaximaCounter2      << '\t' // [11]
       << uniqueMinimaCounter2      << '\t' // [12]
       << period2                   << '\t' // [13]
       << cwIntensity3              << '\t' // [14]
       << uniqueMaximaCounter3      << '\t' // [15]
       << uniqueMinimaCounter3      << '\t' // [16]
       << period3                   << '\t' // [17]
       << deviation01               << '\t' // [18]
       << deviation02               << '\t' // [19]
       << deviation03               << '\t' // [20]
       << deviation12               << '\t' // [21]
       << deviation13               << '\t' // [22]
       << deviation23               << endl;// [23]
  
  cout << np.K_inj                  << '\t' 
       << np.DnuInj                 << '\t' 
       << cwIntensity0              << '\t' 
       << uniqueMaximaCounter0      << '\t' 
       << uniqueMinimaCounter0      << '\t' 
       << period0                   << '\t' 
       << cwIntensity1              << '\t' 
       << uniqueMaximaCounter1      << '\t' 
       << uniqueMinimaCounter1      << '\t' 
       << period1                   << '\t' 
       << cwIntensity2              << '\t' 
       << uniqueMaximaCounter2      << '\t' 
       << uniqueMinimaCounter2      << '\t' 
       << period2                   << '\t' 
       << cwIntensity3              << '\t' 
       << uniqueMaximaCounter3      << '\t' 
       << uniqueMinimaCounter3      << '\t' 
       << period3                   << '\t' 
       << deviation01               << '\t' 
       << deviation02               << '\t' 
       << deviation03               << '\t' 
       << deviation12               << '\t' 
       << deviation13               << '\t' 
       << deviation23               << endl;
  
}



void Analyzer::writeInjectionBifurcationData (ofstream& file)
{
  
  if (fixedPoint0 == true) {
    file << np.DnuInj << '\t' << cwIntensity0 << endl;
    cout << np.DnuInj << '\t' << cwIntensity0 << endl;
  }
  
  else {
    
    for (unsigned int i = 0; i < uniqueMaxima0.size(); i++) {
      file << np.DnuInj << '\t' << uniqueMaxima0[i] << endl;
      cout << np.DnuInj << '\t' << uniqueMaxima0[i] << endl;
    }
    
    for (unsigned int i = 0; i < uniqueMinima0.size(); i++) {
      file << np.DnuInj << '\t' << uniqueMinima0[i] << endl;
      cout << np.DnuInj << '\t' << uniqueMinima0[i] << endl;
    }
    
  }
  
}



void Analyzer::writeFeedbackSingleLaserData (ofstream& file)
{
  
  file << p.K_fb << '\t' << p.phi   << '\t' 
       << cwIntensity0              << '\t' 
       << uniqueMaximaCounter0      << '\t' 
       << uniqueMinimaCounter0      << '\t' 
       << period0                   << '\t' 
       << LyapunovExponent          << endl;
  
  cout << p.K_fb << '\t' << p.phi   << '\t' 
       << cwIntensity0              << '\t' 
       << uniqueMaximaCounter0      << '\t' 
       << uniqueMinimaCounter0      << '\t' 
       << period0                   << '\t' 
       << LyapunovExponent          << endl;
  
}



void Analyzer::writeFeedbackLaserNetworkData (ofstream& file)
{
  
  file << p.K_fb                    << '\t' // [0]
       << p.phi                     << '\t' // [1]
       << cwIntensity0              << '\t' // [2]
       << uniqueMaximaCounter0      << '\t' // [3]
       << uniqueMinimaCounter0      << '\t' // [4]
       << period0                   << '\t' // [5]
       << cwIntensity1              << '\t' // [6]
       << uniqueMaximaCounter1      << '\t' // [7]
       << uniqueMinimaCounter1      << '\t' // [8]
       << period1                   << '\t' // [9]
       << cwIntensity2              << '\t' // [10]
       << uniqueMaximaCounter2      << '\t' // [11]
       << uniqueMinimaCounter2      << '\t' // [12]
       << period2                   << '\t' // [13]
       << cwIntensity3              << '\t' // [14]
       << uniqueMaximaCounter3      << '\t' // [15]
       << uniqueMinimaCounter3      << '\t' // [16]
       << period3                   << '\t' // [17]
       << deviation01               << '\t' // [18]
       << deviation02               << '\t' // [19]
       << deviation03               << '\t' // [20]
       << deviation12               << '\t' // [21]
       << deviation13               << '\t' // [22]
       << deviation23               << endl;// [23]
  
  cout << p.K_fb                    << '\t' 
       << p.phi                     << '\t' 
       << cwIntensity0              << '\t' 
       << uniqueMaximaCounter0      << '\t' 
       << uniqueMinimaCounter0      << '\t' 
       << period0                   << '\t' 
       << cwIntensity1              << '\t' 
       << uniqueMaximaCounter1      << '\t' 
       << uniqueMinimaCounter1      << '\t' 
       << period1                   << '\t' 
       << cwIntensity2              << '\t' 
       << uniqueMaximaCounter2      << '\t' 
       << uniqueMinimaCounter2      << '\t' 
       << period2                   << '\t' 
       << cwIntensity3              << '\t' 
       << uniqueMaximaCounter3      << '\t' 
       << uniqueMinimaCounter3      << '\t' 
       << period3                   << '\t' 
       << deviation01               << '\t' 
       << deviation02               << '\t' 
       << deviation03               << '\t' 
       << deviation12               << '\t' 
       << deviation13               << '\t' 
       << deviation23               << endl;
  
}



void Analyzer::writeFeedbackBifurcationData (ofstream& file)
{
  
  if (fixedPoint0 == true) {
    file << p.K_fb << '\t' << cwIntensity0 << endl;
    cout << p.K_fb << '\t' << cwIntensity0 << endl;
  }
  
  else {
    
    for (unsigned int i = 0; i < uniqueMaxima0.size(); i++) {
      file << p.K_fb << '\t' << uniqueMaxima0[i] << endl;
      cout << p.K_fb << '\t' << uniqueMaxima0[i] << endl;
    }
    
    for (unsigned int i = 0; i < uniqueMinima0.size(); i++) {
      file << p.K_fb << '\t' << uniqueMinima0[i] << endl;
      cout << p.K_fb << '\t' << uniqueMinima0[i] << endl;
    }
    
  }
  
}

