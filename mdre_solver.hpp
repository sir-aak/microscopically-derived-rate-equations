#ifndef MDRE_SOLVER_H
#define MDRE_SOLVER_H


class Solver
{
  
  private:
  
  int    i;
  double tol;
  
  // Mersenne twister random number generator for stochastic integrals
  mt19937 rng;
  
  normal_distribution<double> WienerIncrement;
  
  
  // variables for deterministic methods
  
  arma::cx_vec X_temp;
  arma::cx_vec X_temp_prev;
  
  arma::cx_vec k1;
  arma::cx_vec k2;
  arma::cx_vec k3;
  arma::cx_vec k4;
  
  arma::colvec dk1;
  arma::colvec dk2;
  arma::colvec dk3;
  arma::colvec dk4;
  
  // variables for stochastic methods
  
  arma::cx_vec DW;
  
  arma::cx_vec a;
  arma::cx_vec b;
  
  
  
  public:
  
  Solver ();
  
  arma::cx_vec dW (Network& N);
  
  void explicitEuler (Network& N, double t);
  void implicitEuler (Network& N, double t, double _tol = 1e-10);
  void explicitMidpoint (Network& N, double t);
  void explicitHeun (Network& N, double t);
  void explicitRungeKutta4 (Network& N, double t);
  
  void explicitEulerDelay (Network& N, double t);
  void explicitRungeKutta4Delay (Network& N, double t);
  
  void explicitEulerMaruyama (Network& N, double t);
  void explicitEulerMaruyamaDelay (Network& N, double t);
  
  void explicitEulerLinearized (Network& N, double t);
  void explicitRungeKutta4Linearized (Network& N, double t);
  void explicitEulerDelayLinearized (Network& N, double t);
  void explicitRungeKutta4DelayLinearized (Network& N, double t);
  
};


#endif

