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
};

vector<double> mediarannyu(double numblocchi, double throws, Random rand);
vector<double> mediarannyuminmax(double numblocchi, double throws, Random rand, double min, double max);
vector<double> mediaexp(double numblocchi, double throws, Random rand, double mean);
vector<double> mediaCL(double numblocchi, double throws, Random rand, double gamma, double mu);
vector<double> chi2(double interv, double n, double num_iterations, Random rnd);
void generateData(const string& filename, int throws, int N, Random& rnd);
vector<vector<double>> elaboraDati(double numblocchi, double throws, const vector<double>& media);
void scriviSuFile(string nomeFile, const vector<vector<double>>& vettori);
vector<double> calcolaPi(int nblocchi, int L, double d, double lungh, Random& rnd);
vector <double> calcolaPi_hom(int nblocchi, int L, double d, double lungh, Random& rnd);
vector <double> calcolaPi_con_pi(int nblocchi, int L, double d, double lungh, Random& rnd);


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
