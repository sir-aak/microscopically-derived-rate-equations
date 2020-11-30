#ifndef MDRE_TIME_SERIES_H
#define MDRE_TIME_SERIES_H


class TimeSeries
{
  
  private:
  
  mt19937 rng;
  uniform_real_distribution<double> uniform;
  
  public:
  
  double  t;
  Network N;
  Solver  S;
  
  arma::colvec timeValues;
  arma::colvec intensityValues0;
  arma::colvec intensityValues1;
  arma::colvec intensityValues2;
  arma::colvec intensityValues3;
  arma::colvec LyapunovValues;
  
  vector<string> deterministic_methods;
  vector<string> stochastic_methods;
  vector<string> no_delay_methods;
  vector<string> delay_methods;
  
  arma::cx_vec dmdre;
  
  double drho_GS_e_act;
  double drho_GS_h_act;
  double drho_GS_e_inact;
  double drho_GS_h_inact;
  double drho_ES_e;
  double drho_ES_h;
  double dw_e;
  double dw_h;
  
  
  TimeSeries ();
  TimeSeries (int vectorLength);
  
  void   updateNumericParameters ();
  void   setIntegrationMethods ();
  void   checkTimeSeriesFeasibility ();
  void   setRotatingFrame ();
  double electronSum ();
  double holeSum ();
  double chargeConservation ();
  double dchargeConservation ();
  void   timeSeriesStep (Network& _N, long long int i);
  void   writeTimeSeries (ostream& file);
  void   networkPerturbation ();
  void   showParameters ();
  
};


#endif

