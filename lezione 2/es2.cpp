#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "random.h"

using namespace std;

int main() {
    // Inizializzazione oggetto Random per la generazione di numeri casuali
    Random rnd;

    // Definizione dei parametri della simulazione
    double throws = 10000;
    double passi = 100;
    double blocchi = 100;
    double a = 1;

    // Esecuzione del random walk discreto e continuo
    vector<vector<double>> risultati_discr = random_walk_average(passi, throws, blocchi, a, rnd);
    vector<vector<double>> risultati_cont = random_walk_average_cont(passi, throws, blocchi, a, rnd);

    // Scrittura dei risultati su file
    scriviSuFile("dati.dat", {risultati_discr[0], risultati_discr[2], risultati_cont[0], risultati_cont[2]});

    return 0;
}

