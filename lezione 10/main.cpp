#include <iostream>
#include <vector>
#include <limits>
#include "mpi.h"
#include "tsp.h"
#include "random.h"

using namespace std;

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int numCities = 110;
    Random rand;
    rand.initialize(rank); // Initialize the random number generator with MPI rank
    vector<City> cities = createCitiesFromDataFile("cap_prov_ita.dat");

    GeneticAlgorithm ga(cities, 1000, 0.2, 0.9, rand);
    ga.initializePopulation();

    int gen = 1000;
    int M = 10; // Synchronization frequency

    double bestFitness = numeric_limits<double>::max();
    Individual bestIndividual{vector<int>(numCities)};
    vector<double> fitnessPerNode(size);
    string fitnessFilename = "fitness_per_generation.dat";

    for (int i = 0; i < gen; i++) {
        ga.evolve();
        Individual localBest = ga.getBestIndividual();

        if (localBest.fitness < bestFitness) {
            bestFitness = localBest.fitness;
            bestIndividual = localBest;
        }

        if (i % M == 0) {
            vector<int> localBestPath = bestIndividual.path;
            vector<int> globalBestPath(numCities);
            vector<int> allPaths(size * numCities);

            MPI_Gather(localBestPath.data(), numCities, MPI_INT, allPaths.data(), numCities, MPI_INT, 0, MPI_COMM_WORLD);

            if (rank == 0) {
                double globalBestFitness = bestFitness;
                vector<int> globalBestIndividualPath = localBestPath;

                for (int p = 1; p < size; p++) {
                    vector<int> tmpPath(allPaths.begin() + p * numCities, allPaths.begin() + (p + 1) * numCities);
                    Individual tmpIndividual(tmpPath);
                    tmpIndividual.calculateFitness(cities);

                    if (tmpIndividual.fitness < globalBestFitness) {
                        globalBestFitness = tmpIndividual.fitness;
                        globalBestIndividualPath = tmpPath;
                    }
                }

                bestFitness = globalBestFitness;
                bestIndividual.path = globalBestIndividualPath;
                bestIndividual.calculateFitness(cities);
            }

            MPI_Bcast(bestIndividual.path.data(), numCities, MPI_INT, 0, MPI_COMM_WORLD);
        }

        // Save fitness per generation
        fitnessPerNode[rank] = bestFitness;
        MPI_Gather(&bestFitness, 1, MPI_DOUBLE, fitnessPerNode.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            saveFitnessPerGeneration(fitnessFilename, fitnessPerNode, i);
        }
    }

    if (rank == 0) {
        cout << "Best path found: " << endl;
        for (int cityIndex : bestIndividual.path) {
            cout << cityIndex << " ";
        }
        cout << endl;
        cout << "Best path length: " << bestIndividual.fitness << endl;
    }
    saveBestIndividualToFile("italy.dat", bestIndividual, cities);

    MPI_Finalize();
    return 0;
}
