#ifndef MDRE_AUXILIARY_H
#define MDRE_AUXILIARY_H


#include "mdre_time_series.hpp"


void updateNumericParameters ();

arma::rowvec nonlinearScatteringRatesVector (TimeSeries& TS, double w);

arma::rowvec scatteringProcessesVector (TimeSeries& TS, double w);

double tau_e_inverse (double rho_ES_e, double w_e);

double tau_h_inverse (double rho_ES_h, double w_h);

void createArgumentFile (string argfileType);


#endif

