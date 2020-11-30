#ifndef MDRE_TISI_FINDER_H
#define MDRE_TISI_FINDER_H 


class TISI_Finder
{
  
  public:
  
  int    numberOfIntervals;
  double tdead;
  double tlimit0, tlimit1, tlimit2, tlimit3;
  double spikeThreshold_up;
  double spikeThreshold_low;
  bool   wait0, wait1, wait2, wait3;
  
  // vectors for storage of all times at which spikes occur
  arma::colvec spikeTimes0;
  arma::colvec spikeTimes1;
  arma::colvec spikeTimes2;
  arma::colvec spikeTimes3;
  
  // vectors for storage of inter-spike intervall times
  arma::colvec Tisi0;
  arma::colvec Tisi1;
  arma::colvec Tisi2;
  arma::colvec Tisi3;
  
  
  TISI_Finder ();
  
  TISI_Finder (int _numberOfIntervals);
  
  double roundDown (double a, int decimal);
  
  double roundUp (double a, int decimal);
  
  double findSNIPER (double bifurcationTol = 1e-15, 
                     double intensityTol = 1e-2, bool down = true);
  
  double findHomoclinic (double bifurcationTol = 1e-7, 
                         double intensityTol = 1e-2, bool left = true);
  
  void storeSpikeTimes (int& spikeCounter0, arma::colvec& timeValues, 
                        arma::colvec& intensityValues0);
  
  void storeSpikeTimes (int& spikeCounter0, int& spikeCounter1, 
                       int& spikeCounter2, int& spikeCounter3, 
                       arma::colvec& timeValues, 
                       arma::colvec& intensityValues0, 
                       arma::colvec& intensityValues1, 
                       arma::colvec& intensityValues2, 
                       arma::colvec& intensityValues3);
  
  void getTISI ();
  
  void showSpikeTimesAndTISI ();
  
  void writeCumulantsToFile (ofstream& file, double t);
  
  void writeCumulantsToFileLaserNetwork (ofstream& file, double t);
  
};


#endif

