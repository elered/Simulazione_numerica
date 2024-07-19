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


class funz {

public:

    double eval(int func_index, double x) const {
        switch (func_index) {
            case 1:
                return eval_gx_1(x);
            case 2:
                return eval_px_1(x);
            case 3:
                return eval_gx_2(x);
            case 4:
                return eval_px_2(x);
            case 5:
                return eval_fx(x);
            default:
                return 0.0; // Gestione errore
        }
    }

private:

    double eval_gx_1(double x) const {
        return (M_PI / 3) * ((cos((M_PI * x) / 2)) / (1 - (x * x)));
    }

    double eval_px_1(double x) const {
        return ((3 / M_PI) * (1 - (x * x)));
    }

    double eval_gx_2(double x) const {
        return (M_PI / 2) * cos((M_PI * x) / 2) * ((3 * M_PI - 2) / (6 * (-x * x + M_PI / 2)));
    }

    double eval_px_2(double x) const {
        return (6 * (-x * x + M_PI / 2)) / (3 * M_PI - 2);
    }

    double eval_fx(double x) const {
        return (M_PI / 2) * cos((M_PI * x) / 2);
    }
};


vector<vector<double>> random_walk_average(int num_steps, int num_trials, int num_blocks, double a, Random &rnd);
vector<vector<double>> random_walk_average_cont(int num_steps, int num_trials, int num_blocks, double a, Random &rnd);
vector<double> mediarannyu(double numblocchi, double throws, Random rand);
vector<double> integraleave(double xmin, double xmax, int throws, int blocchi,Random &rnd, funz *f, int funz);
void scriviSuFile(string nomeFile, const vector<vector<double>>& vettori);
vector<vector<double>> elaboraDati(double numblocchi, double throws, const vector<double>& media);
vector<double> integralehom(double xmin, double xmax, double fmax, Random &rnd,int throws, int blocchi, funz *f, int prob, int funzione);

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
