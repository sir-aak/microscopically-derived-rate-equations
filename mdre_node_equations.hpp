#ifndef MDRE_NODE_EQUATIONS_H
#define MDRE_NODE_EQUATIONS_H


class Node
{
  
  private:
  
  arma::cx_vec dX_node;
  arma::cx_vec dX_node_stoch;
  arma::colvec dX_node_transformed;
  
  
  public:
  
  complex<double> E;
  double rho_GS_e_act;
  double rho_GS_h_act;
  double rho_GS_e_inact;
  double rho_GS_h_inact;
  double rho_ES_e;
  double rho_ES_h;
  double w_e;
  double w_h;
  
  double d_Ex;
  double d_Ey;
  double d_rho_GS_e_act;
  double d_rho_GS_h_act;
  double d_rho_GS_e_inact;
  double d_rho_GS_h_inact;
  double d_rho_ES_e;
  double d_rho_ES_h;
  double d_w_e;
  double d_w_h;
  
  int          nodeNumber;
  arma::cx_vec X_node;
  arma::colvec X_node_transformed;
  
  Node (int _nodeNumber, arma::cx_vec& X_net);
  
  double getA ();
  double getphi ();
  arma::cx_vec get_X_node (arma::cx_vec& X_net);
  arma::cx_vec f_node (arma::cx_vec& X_net, double t);
  arma::cx_vec f_node_stoch ();
  
  double g ();
  double delta_omega ();
  double D ();
  complex<double> dE_sp ();
  
  double E_F_e_eq ();
  double E_F_h_eq ();
  
  double R_sp_GS_act ();
  double R_sp_GS_inact ();
  double R_sp_ES ();
  double r_loss_w ();
  
  double S_GS_e_cap_in ();
  double S_GS_h_cap_in ();
  double S_ES_e_cap_in ();
  double S_ES_h_cap_in ();
  double S_e_rel_in ();
  double S_h_rel_in ();
  
  double S_GS_e_cap_out ();
  double S_GS_h_cap_out ();
  double S_ES_e_cap_out ();
  double S_ES_h_cap_out ();
  double S_e_rel_out ();
  double S_h_rel_out ();
  
  double S_GS_e_cap_act ();
  double S_GS_h_cap_act ();
  double S_GS_e_cap_inact ();
  double S_GS_h_cap_inact ();
  double S_ES_e_cap ();
  double S_ES_h_cap ();
  
  double S_e_rel_act ();
  double S_h_rel_act ();
  double S_e_rel_inact ();
  double S_h_rel_inact ();
  
  
  arma::colvec get_X_node_transformed (arma::colvec& X_net_transformed);
  arma::colvec f_node_transformed (arma::colvec& X_net, double t);
  
  double d_g ();
  double d_delta_omega ();
  double d_D ();
  
  double d_E_F_e_eq ();
  double d_E_F_h_eq ();
  
  double d_R_sp_GS_act ();
  double d_R_sp_GS_inact ();
  double d_R_sp_ES ();
  double d_r_loss_w ();
  
  double d_S_GS_e_cap_in ();
  double d_S_GS_h_cap_in ();
  double d_S_ES_e_cap_in ();
  double d_S_ES_h_cap_in ();
  double d_S_e_rel_in ();
  double d_S_h_rel_in ();
  
  double d_S_GS_e_cap_out ();
  double d_S_GS_h_cap_out ();
  double d_S_ES_e_cap_out ();
  double d_S_ES_h_cap_out ();
  double d_S_e_rel_out ();
  double d_S_h_rel_out ();
  
  double d_S_GS_e_cap_act ();
  double d_S_GS_h_cap_act ();
  double d_S_GS_e_cap_inact ();
  double d_S_GS_h_cap_inact ();
  double d_S_ES_e_cap ();
  double d_S_ES_h_cap ();
  
  double d_S_e_rel_act ();
  double d_S_h_rel_act ();
  double d_S_e_rel_inact ();
  double d_S_h_rel_inact ();
  
};


#endif

