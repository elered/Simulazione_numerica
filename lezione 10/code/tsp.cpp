#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>
#include <unordered_set>
#include "tsp.h"
#include "random.h"

using namespace std;

City::City(double x, double y) : x(x), y(y) {}

double City::distanceTo(const City& other_city) const {
    return sqrt(pow(x - other_city.x, 2) + pow(y - other_city.y, 2));
}

Individual::Individual(const vector<int>& path) : path(path), fitness(0) {}

void Individual::calculateFitness(const vector<City>& cities) {
    fitness = 0;
    for (size_t i = 0; i < path.size() - 1; i++) {
        fitness += cities[path[i]].distanceTo(cities[path[i + 1]]);
    }
    fitness += cities[path.back()].distanceTo(cities[path[0]]); // return to start
}

bool Individual::isValid(int numCities) {
    unordered_set<int> uniqueCities(path.begin(), path.end());
    return uniqueCities.size() == static_cast<std::unordered_set<int>::size_type>(numCities);
}

GeneticAlgorithm::GeneticAlgorithm(const vector<City>& cities, size_t populationSize, double mutationRate, double crossoverRate, Random &rnd)
    : cities(cities), populationSize(populationSize), mutationRate(mutationRate), crossoverRate(crossoverRate) {
    initializePopulation();
    rand = rnd;
}

void GeneticAlgorithm::initializePopulation() {
    vector<int> basePath(cities.size());
    iota(basePath.begin(), basePath.end(), 0); // initialize base path with 1, ..., N

    for (int i = 0; i < populationSize; i++) {
        random_shuffle(basePath.begin() + 1, basePath.end()); // Random shuffle
        Individual individual(basePath);
        individual.calculateFitness(cities);
        population.push_back(individual);
    }
}

Individual GeneticAlgorithm::selectParent(double p) {
    sort(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
        return a.fitness > b.fitness; // Assume higher fitness is better
    });

    double r = rand.Rannyu(0, 1);
    int M = population.size();
    int j = static_cast<int>(M * pow(r, p));

    return population[j];
}

void GeneticAlgorithm::shiftCities(Individual& individual, int start, int m, int n) {
    if (rand.Rannyu(0, 1) < mutationRate) {
        int N = individual.path.size();

        if (start < 1 || start >= N || m <= 0 || start + m > N) return;

        vector<int> temp(individual.path.begin() + start, individual.path.begin() + start + m);

        int shiftIndex = (start + n) % N;
        if (shiftIndex < 1) shiftIndex += N;
        if (shiftIndex >= start && shiftIndex < start + m) return;

        individual.path.erase(individual.path.begin() + start, individual.path.begin() + start + m);

        if (shiftIndex > start) {
            shiftIndex -= m; // Adjust shiftIndex to reflect the removal of elements before it
        }

        individual.path.insert(individual.path.begin() + shiftIndex, temp.begin(), temp.end());
    }
}

void GeneticAlgorithm::invertCities(Individual& individual, int start, int m) {
    if (rand.Rannyu(0, 1) < mutationRate) {
        int N = individual.path.size();
        if (start < 1 || start + m > N) return;

        reverse(individual.path.begin() + start, individual.path.begin() + start + m);
    }
}

void GeneticAlgorithm::mutate(Individual& individual) {
    if (rand.Rannyu(0, 1) < mutationRate) {
        int idx1 = static_cast<int>(rand.Rannyu(1, cities.size()));
        int idx2 = static_cast<int>(rand.Rannyu(1, cities.size()));
        while (idx1 == idx2) {
            idx2 = static_cast<int>(rand.Rannyu(1, cities.size()));
        }

        swap(individual.path[idx1], individual.path[idx2]);
    }
}

void GeneticAlgorithm::crossover(const Individual& parent1, const Individual& parent2, Individual& offspring1, Individual& offspring2) {
    if (rand.Rannyu(0, 1) < crossoverRate) {
        int cut = static_cast<int>(rand.Rannyu(1, parent1.path.size() - 1));
        offspring1.path = vector<int>(parent1.path.begin(), parent1.path.begin() + cut);
        offspring2.path = vector<int>(parent2.path.begin(), parent2.path.begin() + cut);

        auto addRemainingCities = [](const vector<int>& parentPath, vector<int>& offspringPath) {
            unordered_set<int> citiesInOffspring(offspringPath.begin(), offspringPath.end());
            for (int city : parentPath) {
                if (citiesInOffspring.find(city) == citiesInOffspring.end()) {
                    offspringPath.push_back(city);
                    citiesInOffspring.insert(city);
                }
            }
        };

        addRemainingCities(parent2.path, offspring1.path);
        addRemainingCities(parent1.path, offspring2.path);
    } else {
        offspring1 = parent1;
        offspring2 = parent2;
    }
}

void GeneticAlgorithm::evolve() {
    vector<Individual> newPopulation;

    while (newPopulation.size() < static_cast<std::vector<Individual>::size_type>(populationSize)) {
        Individual parent1 = selectParent(0.2);
        Individual parent2 = selectParent(0.2);

        Individual offspring1(parent1.path);
        Individual offspring2(parent2.path);

        crossover(parent1, parent2, offspring1, offspring2);

        mutate(offspring1);
        mutate(offspring2);

        int startShift = static_cast<int>(rand.Rannyu(1, offspring1.path.size() - 1));
        int mShift = static_cast<int>(rand.Rannyu(1, offspring1.path.size() - startShift));
        int nShift = static_cast<int>(rand.Rannyu(1, offspring1.path.size()));

        shiftCities(offspring1, startShift, mShift, nShift);

        int startInvert = static_cast<int>(rand.Rannyu(1, offspring2.path.size()));
        int mInvert = static_cast<int>(rand.Rannyu(0, offspring2.path.size() - startInvert));

        invertCities(offspring2, startInvert, mInvert);

        offspring1.calculateFitness(cities);
        offspring2.calculateFitness(cities);

        if (offspring1.isValid(cities.size()) && offspring2.isValid(cities.size())) {
            newPopulation.push_back(offspring1);
            newPopulation.push_back(offspring2);
        }
    }

    population = newPopulation;
}

Individual GeneticAlgorithm::getBestIndividual() {
    return *min_element(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
        return a.fitness < b.fitness;
    });
}

void saveBestIndividualToFile(const std::string& filename, const Individual& best, const std::vector<City>& cities) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    for (int cityIndex : best.path) {
        outFile << cities[cityIndex].x << " " << cities[cityIndex].y << "\n";
    }

    outFile.close();
}

vector<City> createCitiesFromDataFile(const std::string& filename) {
    vector<City> cities;
    ifstream inFile(filename);
    if (!inFile) {
        cerr << "Error opening file: " << filename << endl;
        return cities;
    }

    double x, y;
    while (inFile >> x >> y) {
        cities.push_back(City(x, y));
    }

    inFile.close();
    return cities;
}

void saveFitnessPerGeneration(const std::string& filename, const vector<double>& fitnessPerNode, int generation) {
    std::ofstream outFile;
    outFile.open(filename, std::ios_base::app);
    if (!outFile) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    outFile << generation;
    for (const auto& fitness : fitnessPerNode) {
        outFile << " " << fitness;
    }
    outFile << "\n";
    outFile.close();
}
