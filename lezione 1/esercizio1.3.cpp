#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include <vector>

using namespace std;

int main() {
    Random rnd;

    double d = 10;    // distanza tra le righe
    double lungh = 5; // lunghezza ago
    double throws = 100000; // numero di lanci
    double nblocchi = 100; // numero di blocchi
    double counter = 0;
    double xi, xf = 0;
    double dx;
    double L = throws / nblocchi;

    // Calcolo di pi con tre diversi metodi
    vector<double> pi = calcolaPi(nblocchi, L, d, lungh, rnd);
    vector<double> pi_hom = calcolaPi_hom(nblocchi, L, d, lungh, rnd);
    vector<double> pi_con_pi = calcolaPi_con_pi(nblocchi, L, d, lungh, rnd);

    // Elaborazione dei risultati
    vector<vector<double>> risultati_pi = elaboraDati(nblocchi, throws, pi);
    vector<vector<double>> risultati_pi_hom = elaboraDati(nblocchi, throws, pi_hom);
    vector<vector<double>> risultati_pi_con_pi = elaboraDati(nblocchi, throws, pi_con_pi);

    // Scrittura dei risultati su file
    scriviSuFile("pi.dat", {risultati_pi[0], risultati_pi[2], risultati_pi_hom[0], risultati_pi_hom[2], risultati_pi_con_pi[0], risultati_pi_con_pi[2]});

    return 0;
}




