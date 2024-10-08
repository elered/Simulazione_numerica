/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Random__
#define __Random__
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
using namespace std;

// This class contains functions for generating random numbers using the RANNYU algorithm
class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // Default constructor
  //Random();
  // Destructor
  ~Random();
  // Method to set the seed for the RNG
  void SetRandom(int * , int, int);
  // Method to save the seed to a file
  void SaveSeed();
  // Method to generate a random number in the range [0,1)
  Random();
  //Mio costruttore random
  double Rannyu(void);
  // Method to generate a random number in the range [min,max)
  double Rannyu(double min, double max);
  // Method to generate a random number with a Gaussian distribution
  double Gauss(double mean, double sigma);
  //Method to generate a random number with an exponential ditribution
  double Exp(double mean);
  //Method to generate a random number with a Cauchy-Lorentz distribution
  double CauchyLor (double gamma, double mu);

  double Cos();

  double Sin();
};

vector<double> mediarannyu(double numblocchi, double throws, Random rand);
vector<double> integraleave(double xmin, double xmax, int throws, int blocchi, Random& rnd);
vector<double> integralehom(double xmin, double xmax, double fmax, Random& rnd, int throws, int blocchi);
vector <double> GBM(int asset_price, double k, double r, int T, double mean, double sigma, int numblocchi, int throws, Random rand, int call);
vector <double> GBM_geom (int asset_price,double k, double r, int T, double mean, double sigma, int numblocchi, int throws, Random rand, int call);
vector<vector<double>> elaboraDati(double numblocchi, double throws, const vector<double>& media);
void scriviSuFile(string nomeFile, const vector<vector<double>>& vettori);

#endif // __Random__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
