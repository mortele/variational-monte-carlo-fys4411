#include <iostream>
#include <memory>
#include "system.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"

using namespace std;


int main() {
    // Seed for the random number generator
    int seed = 2023;

    unsigned int numberOfDimensions = 1;
    unsigned int numberOfParticles = 1;
    unsigned int numberOfMetropolisSteps = (unsigned int) 1e6;
    double omega = 1.0; // Oscillator frequency.
    double alpha = 0.5; // Variational parameter.
    double stepLength = 0.1; // Metropolis step length.
    double equilibration = 0.1; // Amount of the total steps used for
                                // equilibration.

    // Construct a shared pointer to a new System
    auto system = std::make_shared<System>(seed);
    system->setHamiltonian(std::make_unique<HarmonicOscillator>(system, omega));
    system->setWaveFunction(std::make_unique<SimpleGaussian>(system, alpha));
    system->setInitialState(std::make_unique<RandomUniform>(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction(equilibration);
    system->setStepLength(stepLength);
    system->runMetropolisSteps(numberOfMetropolisSteps);

    return 0;
}
