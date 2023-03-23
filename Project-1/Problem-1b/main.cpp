#include <iostream>
#include <vector>
#include <memory>

#include "../../include/system.h"
#include "../../include/simplegaussian.h"
#include "../../include/harmonicoscillator.h"
#include "../../include/initialstate.h"
#include "../../include/metropolis.h"
#include "../../include/random.h"
#include "../../include/particle.h"
#include "../../include/sampler.h"

using namespace std;

int main(int argc, char **argv)
{
    // Seed for the random number generator
    int seed = 2023;

    size_t numberOfDimensions = argc > 1 ? stoi(argv[1]) : 1;
    size_t numberOfParticles = argc > 2 ? stoi(argv[2]) : 1;
    size_t numberOfMetropolisSteps = argc > 3 ? stoi(argv[3]) : 1e6;
    size_t numberOfEquilibrationSteps = numberOfMetropolisSteps/10;
    double omega = 1.0;      // Oscillator frequency.
    double alpha = argc > 4 ? stod(argv[4]) : 0.5;      // Variational parameter.
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
