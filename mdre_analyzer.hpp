#ifndef MDRE_ANALYZER_H
#define MDRE_ANALYZER_H


class Analyzer
{
  
  public:
  
  bool fixedPoint0;
  bool fixedPoint1;
  bool fixedPoint2;
  bool fixedPoint3;
  
  double cwIntensity0;
  double cwIntensity1;
  double cwIntensity2;
  double cwIntensity3;
  
  int maximaCounter0;           int minimaCounter0;
  int maximaCounter1;           int minimaCounter1;
  int maximaCounter2;           int minimaCounter2;
  int maximaCounter3;           int minimaCounter3;
  
  int uniqueMaximaCounter0;     int uniqueMinimaCounter0;
  int uniqueMaximaCounter1;     int uniqueMinimaCounter1;
  int uniqueMaximaCounter2;     int uniqueMinimaCounter2;
  int uniqueMaximaCounter3;     int uniqueMinimaCounter3;
  
  vector<double> maxima0;       vector<double> maximaTime0;
  vector<double> maxima1;       vector<double> maximaTime1;
  vector<double> maxima2;       vector<double> maximaTime2;
  vector<double> maxima3;       vector<double> maximaTime3;
  
  vector<double> minima0;       vector<double> minimaTime0;
  vector<double> minima1;       vector<double> minimaTime1;
  vector<double> minima2;       vector<double> minimaTime2;
  vector<double> minima3;       vector<double> minimaTime3;
  
  vector<double> uniqueMaxima0; vector<double> uniqueMinima0;
  vector<double> uniqueMaxima1; vector<double> uniqueMinima1;
  vector<double> uniqueMaxima2; vector<double> uniqueMinima2;
  vector<double> uniqueMaxima3; vector<double> uniqueMinima3;
  
  double period0;
  double period1;
  double period2;
  double period3;
  
  double deviation01;
  double deviation02;
  double deviation03;
  double deviation12;
  double deviation13;
  double deviation23;
  
  double LyapunovExponent;
  
  
  Analyzer();
  
  void resetData ();
  
  void fixedPointCheck (arma::colvec& intensityValues0, 
                        double tol = 1e-3);
  
  void fixedPointCheck (arma::colvec& intensityValues0, 
                        arma::colvec& intensityValues1, 
                        arma::colvec& intensityValues2, 
                        arma::colvec& intensityValues3, 
                        double tol = 1e-3);
  
  void setcwIntensity (arma::colvec& intensityValues0);
  
  void setcwIntensity (arma::colvec& intensityValues0, 
                       arma::colvec& intensityValues1, 
                       arma::colvec& intensityValues2, 
                       arma::colvec& intensityValues3);
  
  bool isMaximum (double I_prev, double I_curr, double I_next, 
                  double tol = 1e-10);
  
  bool isMinimum (double I_prev, double I_curr, double I_next, 
                  double tol = 1e-10);
  
  bool extremumExists (double I, vector<double> uniqueExtrema, 
                       double tol = 1e-2);
  
  double getPeriod (vector<double> maxima, vector<double> maximaTime, 
                    double tol = 1e-2);
  
  void countAndStoreExtrema (int k, arma::colvec& timeValues, 
                             arma::colvec& intensityValues0);
  
  void countAndStoreExtrema (int k, arma::colvec& timeValues, 
       arma::colvec& intensityValues0, arma::colvec& intensityValues1, 
       arma::colvec& intensityValues2, arma::colvec& intensityValues3);
  
  void analyze (arma::colvec& timeValues, arma::colvec& intensityValues0);
  
  void analyze (arma::colvec& timeValues, arma::colvec& intensityValues0, 
                arma::colvec& LyapunovValues);
  
  void analyze (arma::colvec& timeValues, 
       arma::colvec& intensityValues0, arma::colvec& intensityValues1, 
       arma::colvec& intensityValues2, arma::colvec& intensityValues3);
  
  void writeInjectionSingleLaserData (ofstream& file);
  
  void writeInjectionLaserNetworkData (ofstream& file);
  
  void writeInjectionBifurcationData (ofstream& file);
  
  void writeFeedbackSingleLaserData (ofstream& file);
  
  void writeFeedbackLaserNetworkData (ofstream& file);
  
  void writeFeedbackBifurcationData (ofstream& file);
  
};


#endif

