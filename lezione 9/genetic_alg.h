#ifndef GENETIC_ALG_H
#define GENETIC_ALG_H

#include <vector>
#include <unordered_set>
#include "random.h"

class City {
public:
    double x, y;
    City(double x, double y);
    double distanceTo(const City& other_city) const;
};

class Individual {
public:
    std::vector<int> path;
    double fitness;
    Individual(const std::vector<int>& path);
    void calculateFitness(const std::vector<City>& cities);
    bool isValid(int numCities);
};

class GeneticAlgorithm {
public:
    std::vector<City> cities;
    std::vector<Individual> population;
    int populationSize;
    double mutationRate;
    double crossoverRate;
    Random rand;

    GeneticAlgorithm(const std::vector<City>& cities, size_t populationSize, double mutationRate, double crossoverRate);
    void initializePopulation();
    Individual selectParent(double p);
    void shiftCities(Individual& individual, int start, int m, int n);
    void invertCities(Individual& individual, int start, int m);
    void mutate(Individual& individual);
    void crossover(const Individual& parent1, const Individual& parent2, Individual& offspring1, Individual& offspring2);
    void evolve();
    Individual getBestIndividual();
};

std::vector<City> createCitiesOnCircumference(int numCities, double radius, Random &rand);
std::vector<City> createCitiesInSquare(int numCities, double sideLength, Random& rand);
void saveBestFitnessToFile(const std::string& filename, const std::vector<double>& bestFitnessPerGeneration);
void saveBestIndividualToFile(const std::string& filename, const Individual& best, const std::vector<City>& cities);

#endif 
