#ifndef MDRE_FEATURES_H
#define MDRE_FEATURES_H


void timeSeries ();

void inputOutput_mdre ();

void nonlinearScatteringRates ();

void scatteringProcesses ();

void injectionLineScan ();

void injectionBifurcationLineScan ();

void bifurcationScan (string bifurcationType, double DnuInjCrit);

void feedbackLineScan ();

void feedbackBifurcationLineScan ();

void rotatingFrameChange ();

void SNIPER_change (double bifurcationTol = 1e-7, 
                    double intensityTol = 1e-2, bool down = true);

void homoclinic_change (double bifurcationTol  = 1e-7, 
                        double intensityTol = 1e-2, bool left = true);

void noiseCheck (int numberOfSamples = 1e6);

void TisiVarianceInjection (int numberOfIntervals = 1000);

void TisiVarianceFeedback (int numberOfIntervals = 1000);


#endif

