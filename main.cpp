#include <iostream>
#include <vector>
#include <memory>

#include "system.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "Solvers/metropolis.h"
#include "Math/random.h"
#include "particle.h"
#include "sampler.h"

using namespace std;


int main() {
    // Seed for the random number generator
    int seed = 2023;

    unsigned int numberOfDimensions = 1;
    unsigned int numberOfParticles = 1;
    unsigned int numberOfMetropolisSteps = (unsigned int) 1e6;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1e5;
    double omega = 1.0; // Oscillator frequency.
    double alpha = 0.5; // Variational parameter.
    double stepLength = 0.1; // Metropolis step length.

    // The random engine can also be built without a seed
    auto rng = std::make_unique<Random>(seed);
    // Initialize particles
    auto particles = setupRandomUniformInitialState(stepLength, numberOfDimensions, numberOfParticles, *rng);
    // Construct a unique pointer to a new System
    auto system = std::make_unique<System>(
            // Construct unique_ptr to Hamiltonian
            std::make_unique<HarmonicOscillator>(omega),
            // Construct unique_ptr to wave function
            std::make_unique<SimpleGaussian>(alpha),
            // Construct unique_ptr to solver, and move rng
            std::make_unique<Metropolis>(std::move(rng)),
            // Move the vector of particles to system
            std::move(particles));

    // Run steps to equilibrate particles
    auto acceptedEquilibrationSteps = system->runEquilibrationSteps(
            stepLength,
            numberOfEquilibrationSteps);

    // Run the Metropolis algorithm
    auto sampler = system->runMetropolisSteps(
            stepLength,
            numberOfMetropolisSteps);

    // Output information from the simulation
    sampler->printOutputToTerminal(*system);

    return 0;
}
