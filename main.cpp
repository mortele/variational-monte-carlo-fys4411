#include <iostream>
#include <string>

#include <vector>
#include <memory>
#include <cmath>

#include "system.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "Solvers/metropolis.h"
#include "Math/random.h"
#include "particle.h"
#include "sampler.h"

using namespace std;


int main(int argv, char** argc) {
    // Seed for the random number generator
    int seed = 2023;

    unsigned int numberOfDimensions = 3;
    unsigned int numberOfParticles = 10;
    unsigned int numberOfMetropolisSteps = (unsigned int) 1e6;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1e6;
    double omega = 1.0; // Oscillator frequency.
    double alpha = omega/2.0; // Variational parameter.
    double stepLength = 0.1; // Metropolis step length.
    bool analytical = true;
    string filename = "";

    if( argv == 1 ) {
        cout << "Hello! Usage:" <<endl;
        cout << "./vmc #dims #particles #log10(metropolis-steps) #log10(equilibriation-steps) omega alpha stepLength analytical? filename" << endl;
        cout << "#dims, int: Number of dimensions" << endl;
        cout << "#particles, int: Number of particles" << endl;
        cout << "#log10(metropolis steps), int/double: log10 of number of steps, i.e. 6 gives 1e6 steps" << endl;
        cout << "#log10(equilibriation-steps), int/double: log10 of number of equilibriation steps, i.e. 6 gives 1e6 steps" << endl;
        cout << "omega, double: Trap frequency" << endl;
        cout << "alpha, double: WF parameter for simple gaussian. Analytical sol alpha = omega/2" << endl;
        cout << "stepLenght, double: How far should I move a particle at each MC cycle?" << endl;
        cout << "analytical?, bool: If the analytical expression should be used. Defaults to true" <<endl;
        cout << "filename, string: If the results should be dumped to a file, give the file name. If none is given, a simple print is performed." <<endl;
        return 0;
    }

    if(argv >= 2)
        numberOfDimensions = (unsigned int) atoi(argc[1]);
    if(argv >= 3)
        numberOfParticles = (unsigned int) atoi(argc[2]);
    if(argv >= 4)
        numberOfMetropolisSteps = (unsigned int) pow(10, atof(argc[3]));
    if(argv >= 5)
        numberOfEquilibrationSteps = (unsigned int) pow(10, atof(argc[4]));
    if(argv >= 6)
        omega = (double) atof(argc[5]);
    if(argv >= 7)
        alpha = (double) atof(argc[6]);
    if(argv >= 8)
        stepLength = (double) atof(argc[7]);
    if(argv >= 9)
        analytical = (bool)atoi(argc[8]);
    if(argv >= 10)
        filename = argc[9];

    // The random engine can also be built without a seed
    auto rng = std::make_unique<Random>(seed);
    // Initialize particles
    auto particles = setupRandomUniformInitialState(stepLength, omega, numberOfDimensions, numberOfParticles, *rng);
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

    if(!analytical)
        system->setWaveFunction(std::make_unique<SimpleGaussianNumerical>(alpha));

    // Run steps to equilibrate particles
    auto acceptedEquilibrationSteps = system->runEquilibrationSteps(
            stepLength,
            numberOfEquilibrationSteps);

    // Run the Metropolis algorithm
    auto sampler = system->runMetropolisSteps(
            stepLength,
            numberOfMetropolisSteps);

    // Output information from the simulation
    if(filename == "") {
        sampler->printOutputToTerminal(*system);
    }
    else {
        sampler->writeOutToFile(*system, filename, omega);
    }

    return 0;
}
