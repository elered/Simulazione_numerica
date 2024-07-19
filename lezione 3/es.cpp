#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "random.h"

using namespace std;

int main() {
    // Definizione dei parametri del modello di Black-Scholes-Merton
    int S_0 = 100;       // Prezzo dell'asset al tempo t=0
    int T = 1;            // Tempo di scadenza
    double k = 100;       // Prezzo di esercizio
    double r = 0.1;       // Tasso di interesse risk-free
    double sigma = 0.25;  // Volatilità
    double mean = 0;      // Media
    int blocchi = 100;    // Numero di blocchi per la simulazione
    int lanci = 10000;    // Numero di lanci
    Random rnd;           // Oggetto Random per la generazione di numeri casuali

    // Simulazione del modello di Black-Scholes-Merton per opzioni di tipo Call e Put
    vector<double> gbm_call = GBM(S_0, k, r, T, mean, sigma, blocchi, lanci, rnd, 1); // Utilizzo 1 per la Call
    vector<double> gbm_put = GBM(S_0, k, r, T, mean, sigma, blocchi, lanci, rnd, 0);  // Utilizzo 0 per la Put
    vector<double> gbmgeom_call = GBM_geom(S_0, k, r, T, mean, sigma, blocchi, lanci, rnd, 1);
    vector<double> gbmgeom_put = GBM_geom(S_0, k, r, T, mean, sigma, blocchi, lanci, rnd, 0);

    // Elaborazione dei dati ottenuti dalla simulazione
    vector<vector<double>> risultati_gbm_call = elaboraDati(blocchi, lanci, gbm_call);
    vector<vector<double>> risultati_gbm_put = elaboraDati(blocchi, lanci, gbm_put);
    vector<vector<double>> risultati_gbmgeom_call = elaboraDati(blocchi, lanci, gbmgeom_call);
    vector<vector<double>> risultati_gbmgeom_put = elaboraDati(blocchi, lanci, gbmgeom_put);

    // Scrittura dei risultati su file
    scriviSuFile("risultati.dat", {risultati_gbm_call[0], risultati_gbm_call[2], risultati_gbm_put[0], risultati_gbm_put[2], risultati_gbmgeom_call[0],
                                    risultati_gbmgeom_call[2], risultati_gbmgeom_put[0], risultati_gbmgeom_put[2]});

    return 0;
}
