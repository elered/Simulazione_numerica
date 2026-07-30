#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "random.h"

using namespace std;

int main() {
    // Inizializzazione dell'oggetto Random
    Random rnd;

    // Definizione dei parametri
    double max = 1;
    double min = 0;
    double throws = 10000;
    double blocchi = 100;

    // Creazione di un'istanza della classe funz tramite puntatore
    funz *funzione = new funz();

    // Calcolo dell'integrale con il metodo della media
    vector<double> ave = integraleave(min, max, throws, blocchi, rnd, funzione, 5);
    vector<vector<double>> risultati_ave = elaboraDati(blocchi, throws, ave);

    // Calcolo dell'integrale con importance sampling usando il metodo eval_px_1
    vector<double> avehom = integralehom(min, max, 3 / 2, rnd, throws, blocchi, funzione, 2, 1);
    vector<vector<double>> risultati_avehom = elaboraDati(blocchi, throws, avehom);

    // Calcolo dell'integrale con importance sampling usando il metodo eval_px_2
    vector<double> avehom_brutto = integralehom(min, max, (3 * M_PI) / (3 * M_PI - 2), rnd, throws, blocchi, funzione, 4, 3);
    vector<vector<double>> risultati_avehom_brutto = elaboraDati(blocchi, throws, avehom_brutto);

    // Scrittura dei risultati su file
    scriviSuFile("risultati.dat", {risultati_ave[0], risultati_ave[2],risultati_avehom[0], risultati_avehom[2],risultati_avehom_brutto[0], risultati_avehom_brutto[2]});

    // Deallocazione della memoria occupata dall'oggetto funzione
    delete funzione;

    return 0;
}
