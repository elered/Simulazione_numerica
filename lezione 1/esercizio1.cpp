#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include <vector>

using namespace std;

int main() {
    Random rnd;

    int M = 100000; // Numero di lanci
    int nblocchi = 100; // Numero di blocchi

    // Chiamata alla funzione media e inizializzazione del vettore Media
    vector<double> Media = mediarannyu(nblocchi, M, rnd);

    vector<vector<double>> risultati_media = elaboraDati(nblocchi, M, Media);

    double L = M / nblocchi; // Numero di lanci in ogni blocco, assicurati che M sia un multiplo di nblocchi
    double num = 0; //Variabile d'appoggio

    // Vettore per le medie dei blocchi 
    vector<double> ave(nblocchi);

    // Calcola le medie dei blocchi
    for (int k = 0; k < nblocchi; k++) {
        num = 0;

        for (int i = 0; i < L; i++) {
            // Calcola il quadrato della differenza tra il numero generato casualmente e 0.5 e sommalo
            num += pow(rnd.Rannyu() - 0.5, 2);
        }   
        // Calcola la media del blocco k
        ave[k] = num / L;
    }

    vector<vector<double>> risultati_ave = elaboraDati(nblocchi, M, ave);

    //Test chi2
    vector<double> vec_chi2 = chi2(100, 10000, 100, rnd);

    // Scrivi i risultati su file
    scriviSuFile("ris.dat", {risultati_media[0], risultati_media[2], risultati_ave[0], risultati_ave[2], vec_chi2});

    return 0;
}

   
