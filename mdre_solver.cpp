#include "mdre_parameters.hpp"
#include "mdre_node_equations.hpp"
#include "mdre_network.hpp"
#include "mdre_solver.hpp"


// link to instantiation of NodeParameters and Parameters class 
// in mdre_processes.cpp
extern NodeParameters np;
extern Parameters p;


// Solver class constructor
Solver::Solver ()
{
  rng.seed(1);
  WienerIncrement = normal_distribution<double>(0.0, p.sqrt_dt);
};



// returns complex vector of Wiener increments
arma::cx_vec Solver::dW (Network& N)
{
  
  DW.zeros(N.X.n_elem);
  
  // single laser case
  if (N.nNodes == 1) {
    DW(0) = WienerIncrement(rng) + 1i * WienerIncrement(rng);
  }
  
  // network case
  else {
    for (int i = 0; i < N.nNodes; i++) {
      DW(9 * i) = WienerIncrement(rng) + 1i * WienerIncrement(rng);
    }
  }
  
  return (DW);
  
}



// methods for deterministic numerical integration of odes:

// performs step with explicit Euler method
void Solver::explicitEuler (Network& N, double t)
{
  N.X = N.X + p.dt * N.f(N.X, t);
}


// performs step with implicit Euler method with fixed point iteration
void Solver::implicitEuler (Network& N, double t, double _tol)
{
  
  tol = _tol;
  
  X_temp      = N.X;
  X_temp_prev = X_temp + 2.0 * tol;
  
  i = 0;
  
  while ((arma::norm(X_temp - X_temp_prev) >= tol) && (i < 15)) {
    
    X_temp_prev = X_temp;
    X_temp      = N.X + p.dt * N.f(X_temp, t);
    i++;
    
  }
  
  N.X = X_temp;
  
}


// performs step with explicit 2-stage / 2nd order Runge-Kutta
// a.k.a. midpoint method
void Solver::explicitMidpoint (Network& N, double t)
{
  
  k1  = N.f(N.X, t);
  k2  = N.f(N.X + p.dt * k1 / 2.0, t + p.dt / 2.0);
  
  N.X = N.X + p.dt * k2;
  
}


// performs step with explicit 2-stage / 2nd order Runge-Kutta
// a.k.a. Heun method
void Solver::explicitHeun (Network& N, double t)
{
  
  k1  = N.f(N.X, t);
  k2  = N.f(N.X + p.dt * k1, t + p.dt);
  
  N.X = N.X + p.dt * (k1 + k2) / 2.0;
  
}


// performs step with explicit 4-stage / 4th order Runge-Kutta method
void Solver::explicitRungeKutta4 (Network& N, double t)
{
  
  k1  = N.f(N.X, t);
  k2  = N.f(N.X + p.dt * k1 / 2.0, t + p.dt / 2.0);
  k3  = N.f(N.X + p.dt * k2 / 2.0, t + p.dt / 2.0);
  k4  = N.f(N.X + p.dt * k3, t + p.dt);
  
  N.X = N.X + p.dt * (k1 + 2.0 * (k2 + k3) + k4) / 6.0;
  
}



// methods for deterministic numerical integration of ddes:

// performs explicit Euler step with delay
void Solver::explicitEulerDelay (Network& N, double t)
{
  
  N.X  = N.X + p.dt * N.f(N.X, N.p0, t);
  N.p0 = N.XHist.col(N.delayedStateIndex);
  N.XHist.col(N.currentStateIndex) = N.X;
  
}


// performs Runge-Kutta 4 step with delay
void Solver::explicitRungeKutta4Delay (Network& N, double t)
{
  
  k1    = N.f(N.X, N.p0, t);
  k2    = N.f(N.X + p.dt * k1 / 2.0, N.p05, t + p.dt / 2.0);
  k3    = N.f(N.X + p.dt * k2 / 2.0, N.p05, t + p.dt / 2.0);
  k4    = N.f(N.X + p.dt * k3, N.p1, t + p.dt);
  
  N.X   = N.X + p.dt * (k1 + 2.0 * (k2 + k3) + k4) / 6.0;
  
  N.p0  = N.XHist.col(N.delayedStateIndex);
  N.p1  = N.XHist.col(N.interpolationIndex);
  N.m0  = N.dXHist.col(N.delayedStateIndex);
  N.m1  = N.dXHist.col(N.interpolationIndex);
  N.p05 = N.Hermite05();
  
  N.XHist.col(N.currentStateIndex)  = N.X;
  N.dXHist.col(N.currentStateIndex) = N.f(N.X, N.p0, t + p.dt);
  
}



// methods for stochastic numerical integration of sdes:

// performs explicit Euler-Maruyama step
// strong convergence order 1/2, weak convergence order 1
void Solver::explicitEulerMaruyama (Network& N, double t)
{
  
  a   = N.f(N.X, t);
  b   = N.f_stoch();
  
  N.X = N.X + a * p.dt + b % dW(N);
  
}



// methods for stochastic numerical integration of sddes

// performs explicit Euler-Maruyama step with delay
void Solver::explicitEulerMaruyamaDelay (Network& N, double t)
{
  
  a    = N.f(N.X, N.p0, t);
  b    = N.f_stoch();
  
  N.X  = N.X + a * p.dt + b % dW(N);
  N.XHist.col(N.currentStateIndex) = N.X;
  N.p0 = N.XHist.col(N.delayedStateIndex);
  
}



// methods for numerical integration of linearized system

// performs step with explicit Euler method for linearized system
void Solver::explicitEulerLinearized (Network& N, double t)
{
  N.deltaX /= norm(N.deltaX);
  N.deltaX  = N.deltaX + p.dt * N.f_linearized(N.deltaX, t);
}


// performs step with explicit 4-stage / 4th order Runge-Kutta method
// for linearized system
void Solver::explicitRungeKutta4Linearized (Network& N, double t)
{
  
  N.deltaX /= norm(N.deltaX);
  
  dk1 = N.f_linearized(N.deltaX, t);
  dk2 = N.f_linearized(N.deltaX + p.dt * dk1 / 2.0, t + p.dt / 2.0);
  dk3 = N.f_linearized(N.deltaX + p.dt * dk2 / 2.0, t + p.dt / 2.0);
  dk4 = N.f_linearized(N.deltaX + p.dt * dk3, t + p.dt);
  
  N.deltaX = N.deltaX + p.dt * (dk1 + 2.0 * (dk2 + dk3) + dk4) / 6.0;
  
}


// performs explicit Euler step with delay for linearized system
void Solver::explicitEulerDelayLinearized (Network& N, double t)
{
  N.deltaX /= norm(N.deltaX);
  N.deltaX  = N.deltaX + p.dt * N.f_linearized(N.deltaX, N.deltap0, t);
  N.deltap0 = N.deltaXHist.col(N.delayedStateIndex);
  N.deltaXHist.col(N.currentStateIndex) = N.deltaX;
}


// performs Runge-Kutta 4 step with delay for linearized system
void Solver::explicitRungeKutta4DelayLinearized (Network& N, double t)
{
  
  N.deltaX /= norm(N.deltaX);
  
  dk1 = N.f_linearized(N.deltaX, N.deltap0, t);
  dk2 = N.f_linearized(N.deltaX + p.dt * dk1 / 2.0, N.deltap05, t + p.dt / 2.0);
  dk3 = N.f_linearized(N.deltaX + p.dt * dk2 / 2.0, N.deltap05, t + p.dt / 2.0);
  dk4 = N.f_linearized(N.deltaX + p.dt * dk3, N.deltap1, t + p.dt);
  
  N.deltaX   = N.deltaX + p.dt * (dk1 + 2.0 * (dk2 + dk3) + dk4) / 6.0;
  
  N.deltap0  = N.deltaXHist.col(N.delayedStateIndex);
  N.deltap1  = N.deltaXHist.col(N.interpolationIndex);
  N.deltam0  = N.ddeltaXHist.col(N.delayedStateIndex);
  N.deltam1  = N.ddeltaXHist.col(N.interpolationIndex);
  N.deltap05 = N.Hermite05Linearized();
  
  N.deltaXHist.col(N.currentStateIndex)  = N.deltaX;
  N.ddeltaXHist.col(N.currentStateIndex) = N.f_linearized(N.deltaX, 
                                             N.deltap0, t + p.dt);
  
}

