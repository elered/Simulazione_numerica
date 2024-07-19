#include <iostream>
#include <vector>
#include <limits>
#include "genetic_alg.h"
#include "random.h"

using namespace std;

int main() {
    int numcities = 34;
    double radius = 1;
    Random rand;
    vector<City> cities = createCitiesOnCircumference(numcities, radius, rand);

    GeneticAlgorithm ga(cities, 1000, 0.2, 0.9);
    ga.initializePopulation();
    double bestFitness = numeric_limits<double>::max();
    Individual bestIndividual{vector<int>(numcities)};
    vector<double> bestFitnessPerGeneration;

    for (int generation = 0; generation < 1000; generation++) {
        ga.evolve();
        Individual best = ga.getBestIndividual();
        bestFitnessPerGeneration.push_back(best.fitness);

        if (best.fitness <= bestFitness) {
            bestFitness = best.fitness;
            bestIndividual = best;
            cout << "Best path length = " << best.fitness << endl;
        }
    }
    saveBestIndividualToFile("best_individual_circle.dat", bestIndividual, cities);
    saveBestFitnessToFile("best_fitness_circle.dat", bestFitnessPerGeneration);

    vector<City> cities_2 = createCitiesInSquare(numcities, 2, rand);
    GeneticAlgorithm ga_2(cities_2, 1000, 0.2, 0.9);
    ga_2.initializePopulation();
    double bestFitness_2 = numeric_limits<double>::max();
    Individual bestIndividual_2{vector<int>(numcities)};
    vector<double> bestFitnessPerGeneration_2;

    for (int generation = 0; generation < 1000; generation++) {
        ga_2.evolve();
        Individual best_2 = ga_2.getBestIndividual();
        bestFitnessPerGeneration_2.push_back(best_2.fitness);

        if (best_2.fitness <= bestFitness_2) {
            bestFitness_2 = best_2.fitness;
            bestIndividual_2 = best_2;
            cout << "Best path length square = " << best_2.fitness << endl;
        }
    }
    saveBestIndividualToFile("best_individual_square.dat", bestIndividual_2, cities_2);
    saveBestFitnessToFile("best_fitness_square.dat", bestFitnessPerGeneration_2);

    return 0;
}
