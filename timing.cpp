#include <iostream>
#include <string>

#include <vector>
#include <memory>
#include <cmath>

#include <chrono>

#include "system.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "Solvers/metropolis.h"
#include "Solvers/metropolishastings.h"
#include "Math/random.h"
#include "particle.h"
#include "sampler.h"

int main(int argv, char **argc)
{
    // Seed for the random number generator
    int seed = 2023;

    // Set default paramters
    unsigned int numberOfDimensions = 3;
    unsigned int numberOfParticles = 10;
    unsigned int numberOfMetropolisSteps = (unsigned int)1e6;
    unsigned int numberOfEquilibrationSteps = (unsigned int)1e6;
    double omega = 1.0;         // Oscillator frequency.
    double alpha = omega / 2.0; // Variational parameter. If using gradient descent, this is the initial guess.
    double stepLength = 0.1;    // Metropolis step length.
    int epochs = 5;             // Number of epochs for gradient descent.
    double lr = 0.1;            // Learning rate for gradient descent.
    double dx = 10e-6;
    bool importanceSampling = false;
    bool gradientDescent = true;
    bool analytical = true;
    double D = 0.5;
    bool timing=0;
    string filename = "";
    
    // If no arguments are given, show usage.
    if (argv == 1)
    {
        cout << "Hello! Usage:" << endl;
        cout << "./vmc #dims #particles #log10(metropolis-steps) #log10(equilibriation-steps) omega alpha stepLength importanceSampling? analytical? gradientDescent? filename" << endl;
        cout << "#dims, int: Number of dimensions" << endl;
        cout << "#particles, int: Number of particles" << endl;
        cout << "#log10(metropolis steps), int/double: log10 of number of steps, i.e. 6 gives 1e6 steps" << endl;
        cout << "#log10(@-steps), int/double: log10 of number of equilibriation steps, i.e. 6 gives 1e6 steps" << endl;
        cout << "omega, double: Trap frequency" << endl;
        cout << "alpha, double: WF parameter for simple gaussian. Analytical sol alpha = omega/2" << endl;
        cout << "stepLenght, double: How far should I move a particle at each MC cycle?" << endl;


        cout << "analytical?, bool: If the analytical expression should be used. Defaults to true" << endl;

        cout << "filename, string: If the results should be dumped to a file, give the file name. If none is given, a simple print is performed." << endl;
            if (argv >= 2)
        numberOfDimensions = (unsigned int)atoi(argc[1]);
    if (argv >= 3)
        numberOfParticles = (unsigned int)atoi(argc[2]);
    if (argv >= 4)
        numberOfMetropolisSteps = (unsigned int)pow(10, atof(argc[3]));
    if (argv >= 5)
        numberOfEquilibrationSteps = (unsigned int)pow(10, atof(argc[4]));
    if (argv >= 6)
        omega = (double)atof(argc[5]);
    if (argv >= 7)
        alpha = (double)atof(argc[6]);
    if (argv >= 8)
        stepLength = (double)atof(argc[7]);
    if (argv >= 9)
        analytical = (bool)atoi(argc[8]);
    if (argv >= 10)
        timing=(bool)atoi(argc[9]);
    if (argv >= 11)
        filename = argc[10];
        // The random engine can also be built without a seed
    auto rng = std::make_unique<Random>(seed);

    // Initialize particles
    auto particles = setupRandomUniformInitialState(stepLength, omega, numberOfDimensions, numberOfParticles, *rng);

    // Construct a unique pointer to a new System
    auto hamiltonian = std::make_unique<HarmonicOscillator>(omega);

    // Initialise SimpleGaussian by default
    std::unique_ptr<class WaveFunction> wavefunction = std::make_unique<SimpleGaussian>(alpha); // Empty wavefunction pointer, since it uses "alpha" in its constructor (can only be moved once).

    // Empty solver pointer, since it uses "rng" in its constructor (can only be moved once).
    std::unique_ptr<class MonteCarlo> solver;

    // Check if numerical gaussian should be used.
    if (!analytical)
        wavefunction = std::make_unique<SimpleGaussianNumerical>(alpha, dx);
    
    solver = std::make_unique<Metropolis>(std::move(rng));
        // Create system pointer, passing in all classes.
    auto system = std::make_unique<System>(
        // Construct unique_ptr to Hamiltonian
        std::move(hamiltonian),
        // Construct unique_ptr to wave function
        std::move(wavefunction),
        // Construct unique_ptr to solver, and move rng
        std::move(solver),
        // Move the vector of particles to system
        std::move(particles));

    // Run steps to equilibrate particles
    auto beginning = std::chrono::high_resolution_clock::now(); //Start timer
    auto acceptedEquilibrationSteps = system->runEquilibrationSteps(
        stepLength,
        numberOfEquilibrationSteps);

    // Run the Metropolis algorithm
    auto sampler = system->runMetropolisSteps(
        stepLength,
        numberOfMetropolisSteps);
    auto ending = std::chrono::high_resolution_clock::now(); //end timer
    double timelapse= std::chrono::duration_cast<std::chrono::microseconds>(ending - beginning).count(); //find time
    // Output information from the simulation, either as file or print
    std::cout<<"timing"<<timing<<std::endl;
    if (filename=="")
    {
        std:cout<<"You need a filename"<<std::endl;
    }
    
    sampler->WriteTimingToFiles(*system, filename, analytical, numberOfEquilibrationSteps, timelapse);

        
        return 0;
    }




}