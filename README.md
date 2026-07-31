# Numerical Simulation Laboratory (NSL)

This repository contains Jupyter Notebooks and Python/C++ scripts developed for the Numerical Simulation Laboratory (NSL) course. The exercises tackle various problems in statistical mechanics, computational physics, and quantitative finance, progressively exploring more advanced numerical techniques.

---

## Exercise 1: Statistical Tests and Monte Carlo Methods

*   **RNG Testing**: Calculation of the mean value and variance for a uniform distribution in the range [0, 1] using the data blocking method.
*   **Chi-Square Test**: Implementation of the $\chi^2$ test to verify the statistical hypothesis that the generated numbers are drawn from a uniform distribution.
*   **Central Limit Theorem**: Verification of the theorem by sampling and adding variables from uniform, exponential, and Cauchy-Lorentz distributions.
*   **Buffon's Experiment**: Stochastic simulation of a needle thrown onto a plane with parallel lines to obtain a numerical estimate of $\pi$.

---

## Exercise 2: Monte Carlo Integration and Random Walk

*   **Monte Carlo Integration**: Evaluation of a 1D integral by sampling a uniform distribution and subsequently applying the Importance Sampling technique to minimize statistical uncertainty.
*   **Random Walk (RW)**: Simulation of a 3D random walk on both a discrete cubic lattice and a continuous space, analyzing the root-mean-square distance from the origin as a function of the number of steps to observe diffusive behavior.

---

## Exercise 3: Plain Vanilla Option Pricing

*   **European Options**: Evaluation of the initial pricing for "Call" and "Put" financial options, comparing the numerical results with the Black-Scholes analytical solution.
*   **Geometric Brownian Motion (GBM)**: Simulation of asset price evolution over time using a GBM, comparing the efficiency of direct sampling of the final price versus a step-by-step discretized sampling.

---

## Exercise 4: Molecular Dynamics

*   **NVE Ensemble Simulation**: Study of the time evolution of an isolated system (constant Number of particles, Volume, and Energy) interacting via a Lennard-Jones potential.
*   **Verlet Algorithm**: Integration of the equations of motion using the Verlet algorithm and application of Periodic Boundary Conditions (PBC) with the minimum image convention to simulate a bulk system.
*   **Phases of Matter**: Equilibration of the thermodynamic system and estimation of fundamental quantities (kinetic, potential, and total energy, temperature, and pressure) for solid, liquid, and gas phases.

---

## Exercise 6: 1D Ising Model

*   **Spin Sampling**: Simulation of a 1D magnetic material based on the Ising model. The system is updated using both the Metropolis algorithm and the Gibbs sampling technique.
*   **Thermodynamics**: Application of the data blocking method to calculate the main macroscopic properties as a function of temperature: internal energy, heat capacity, magnetization (with an external field), and magnetic susceptibility.

---

## Exercise 7: Monte Carlo NVT Simulation

*   **Canonical Ensemble**: Simulation of a Lennard-Jones system in the NVT ensemble (constant Number of particles, Volume, and Temperature) using the Metropolis algorithm.
*   **Autocorrelation and Data Blocking**: Analysis of the autocorrelation function of potential energy and pressure to determine the optimal block size for data blocking.
*   **Radial Distribution Function**: Calculation of the radial distribution function $g(r)$ and comparison of the results with those obtained from the Molecular Dynamics (NVE) code from Exercise 4.

---

## Exercise 8: Variational Monte Carlo (VMC)

*   **1D Quantum Particle**: Application of the Variational Monte Carlo method to find the ground state energy of a single particle in a 1D potential well.
*   **Simulated Annealing**: Implementation of a Simulated Annealing algorithm to optimize the variational parameters of the trial wavefunction and minimize the energy expectation value.

---

## Exercise 9: Genetic Algorithms

*   **Traveling Salesman Problem (TSP)**: Development of a Genetic Algorithm to solve the TSP, finding the shortest path that connects a given set of cities on a plane.
*   **Operators**: Implementation of crossover and various mutation operators (swap, shift, inversion) along with a selection operator based on fitness to evolve the population of paths over generations.

---

## Exercise 10: Simulated Annealing and Parallel Computing

*   **TSP with Simulated Annealing**: Solving the Traveling Salesman Problem using the Simulated Annealing technique, progressively lowering the "temperature" to find the global minimum of the loss function (path length).
*   **MPI Parallelization**: Introduction to parallel computing using MPI (Message Passing Interface) libraries. Running independent searches on different cores (Continents) with periodic migrations of the best paths among nodes to improve the optimization.

---

## Exercise 11: Machine Learning with Neural Networks

*   **Supervised Learning**: Introduction to Machine Learning techniques using Keras and TensorFlow.
*   **Linear and Polynomial Regression**: Training Deep Neural Networks (DNN) to fit linear, polynomial, and trigonometric functions with noisy data.
*   **2D Function Fitting**: Extending the DNN architecture to fit a 2D function, exploring the effects of different activation functions, optimizers, and network topologies (number of layers and neurons).

---

## Exercise 12: Image Recognition with Deep and Convolutional Neural Networks

*   **MNIST Classification with DNN**: Building a Deep Neural Network using Keras to classify handwritten digits (from $0$ to $9$) using the MNIST dataset, which contains $28\times 28$ pixel images[cite: 36]. The data is reshaped, rescaled to a [0, 1] interval, and categorical crossentropy is used alongside the SGD optimizer to train the model over several epochs[cite: 36].
*   **Convolutional Neural Networks (CNN)**: Upgrading the network architecture to a Convolutional Neural Network (CNN)[cite: 36]. This approach utilizes `Conv2D`, `MaxPooling2D`, and `Flatten` layers to take advantage of local spatial correlations and translational invariance in the images, thereby improving classification accuracy[cite: 36].
*   **Custom Image Testing**: Evaluating the trained CNN's real-world performance by testing it on custom handwritten digits created using the GIMP application[cite: 36].
