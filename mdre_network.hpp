#ifndef MDRE_NETWORK_H
#define MDRE_NETWORK_H

#include "mdre_node_equations.hpp"


class Network
{
  
  private:
  
  mt19937 rng;
  uniform_real_distribution<double> uniform;
  arma::cx_vec dX;
  arma::cx_vec dX_stoch;
  
  
  public:
  
  int nNodes;
  
  arma::sp_mat A;
  arma::sp_mat H;
  arma::sp_mat AH;
  arma::sp_mat Ht;
  
  arma::cx_vec ew;
  arma::cx_mat ev;
  
  arma::cx_vec X;
  arma::cx_vec p0;
  arma::cx_vec p1;
  arma::cx_vec m0;
  arma::cx_vec m1;
  arma::cx_vec p05;
  
  arma::cx_mat XHist;
  arma::cx_mat dXHist;
  
  int currentStateIndex;
  int interpolationIndex;
  int delayedStateIndex;
  
  
  arma::colvec deltaX;
  arma::colvec deltap0;
  arma::colvec deltap1;
  arma::colvec deltam0;
  arma::colvec deltam1;
  arma::colvec deltap05;
  
  arma::mat deltaXHist;
  arma::mat ddeltaXHist;
  
  
  vector<Node> nodes;
  
  double dAmp;
  double dphi;
  
  double Amp;
  double Ampd;
  double phi;
  double phid;
  
  arma::colvec Xt;
  arma::mat Eye;
  arma::mat Df_analytic;
  arma::mat Df_numeric;
  
  arma::colvec XEx;
  arma::colvec XEy;
  arma::colvec Xrho_GS_e_act;
  arma::colvec Xrho_GS_h_act;
  arma::colvec Xrho_GS_e_inact;
  arma::colvec Xrho_GS_h_inact;
  arma::colvec Xrho_ES_e;
  arma::colvec Xrho_ES_h;
  arma::colvec Xw_e;
  arma::colvec Xw_h;
  
  arma::colvec fX;
  
  arma::colvec fXEx;
  arma::colvec fXEy;
  arma::colvec fXrho_GS_e_act;
  arma::colvec fXrho_GS_h_act;
  arma::colvec fXrho_GS_e_inact;
  arma::colvec fXrho_GS_h_inact;
  arma::colvec fXrho_ES_e;
  arma::colvec fXrho_ES_h;
  arma::colvec fXw_e;
  arma::colvec fXw_h;
  
  
  arma::colvec XEx_f;
  arma::colvec XEy_f;
  arma::colvec Xrho_GS_e_act_f;
  arma::colvec Xrho_GS_h_act_f;
  arma::colvec Xrho_GS_e_inact_f;
  arma::colvec Xrho_GS_h_inact_f;
  arma::colvec Xrho_ES_e_f;
  arma::colvec Xrho_ES_h_f;
  arma::colvec Xw_e_f;
  arma::colvec Xw_h_f;
  
  arma::colvec XEx_b;
  arma::colvec XEy_b;
  arma::colvec Xrho_GS_e_act_b;
  arma::colvec Xrho_GS_h_act_b;
  arma::colvec Xrho_GS_e_inact_b;
  arma::colvec Xrho_GS_h_inact_b;
  arma::colvec Xrho_ES_e_b;
  arma::colvec Xrho_ES_h_b;
  arma::colvec Xw_e_b;
  arma::colvec Xw_h_b;
  
  arma::colvec fX_f;
  arma::colvec fXEx_f;
  arma::colvec fXEy_f;
  arma::colvec fXrho_GS_e_act_f;
  arma::colvec fXrho_GS_h_act_f;
  arma::colvec fXrho_GS_e_inact_f;
  arma::colvec fXrho_GS_h_inact_f;
  arma::colvec fXrho_ES_e_f;
  arma::colvec fXrho_ES_h_f;
  arma::colvec fXw_e_f;
  arma::colvec fXw_h_f;
  
  arma::colvec fX_b;
  arma::colvec fXEx_b;
  arma::colvec fXEy_b;
  arma::colvec fXrho_GS_e_act_b;
  arma::colvec fXrho_GS_h_act_b;
  arma::colvec fXrho_GS_e_inact_b;
  arma::colvec fXrho_GS_h_inact_b;
  arma::colvec fXrho_ES_e_b;
  arma::colvec fXrho_ES_h_b;
  arma::colvec fXw_e_b;
  arma::colvec fXw_h_b;
  
  
  Network ();
  
  void checkNetworkFeasibility ();
  void setNetworkProperties ();
  void setInitialConditions ();
  void showNetworkEigenSystem ();
  
  arma::cx_vec Hermite05 ();
  
  vector<double> dAmplitude (double t);
  vector<double> dphase (double t);
  
  arma::cx_vec f_local (arma::cx_vec& _X, double t);
  arma::cx_vec f_global (arma::cx_vec& _X);
  arma::cx_vec f_stoch ();
  arma::cx_vec f (arma::cx_vec _X, double t);
  arma::cx_vec f (arma::cx_vec _X, arma::cx_vec _Xd, double t);
  
  
  void setInitialConditionsLinearized ();
  arma::colvec Hermite05Linearized ();
  void setX_transformed ();
  arma::mat JacobianAnalytic (double t);
  arma::mat JacobianNumeric (double t, string differenceType = "central", 
                             double h = 1e-4);
  arma::colvec f_local_linearized (arma::colvec& _deltaX, double t);
  
  arma::colvec f_global_linearized (arma::colvec& _deltaX);
  arma::colvec f_linearized (arma::colvec _deltaX, double t);
  arma::colvec f_linearized (arma::colvec _deltaX, 
                             arma::colvec _deltaXd, double t);
  double LyapunovExponent ();
  
};


#endif

