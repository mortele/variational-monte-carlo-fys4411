#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include "Solvers/metropolis.h"
#include "Solvers/metropolishastings.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/interactinggaussian.h"
#include "particle.h"
#include "sampler.h"
#include "system.h"

using namespace std;

int main(int argv, char **argc)
{
    // Seed for the random number generator
    int seed = 2023;

    // Set default paramters
    unsigned int numberOfDimensions = 3;
    unsigned int numberOfParticles = 10;
    unsigned int numberOfMetropolisSteps = (unsigned int)pow(2, 20);
    unsigned int numberOfEquilibrationSteps = (unsigned int)pow(2, 20);
    double omega = 1.0;                                 // Oscillator frequency.
    double alpha = omega / 2.0;                         // Variational parameter. If using gradient
                                                        // descent, this is the initial guess.
    double stepLength = 0.1;                            // Metropolis step length.
    double epsilon = 0.01;                              // Tolerance for gradient descent.
    double lr = 0.01;                                   // Learning rate for gradient descent.
    double interactionTerm = 0.0043 * sqrt(1. / omega); // Interection constant a. Notice that hbar = m = 1.
    double dx = 10e-6;
    bool importanceSampling = false;
    bool analytical = true;
    bool interacting = false;
    double D = 0.5;
    string filename = "";

    // If no arguments are given, show usage.
    if (argv == 1)
    {
        cout << "Hello! Usage:" << endl;
        cout << "./vmc #dims #particles #log2(metropolis-steps) "
                "#log2(equilibriation-steps) omega alpha stepLength "
                "importanceSampling? analytical? gradientDescent? filename"
             << endl;
        cout << "#dims, int: Number of dimensions" << endl;
        cout << "#particles, int: Number of particles" << endl;
        cout << "#log2(metropolis steps), int/double: log2 of number of steps, i.e. 6 gives 2^6 steps" << endl;
        cout << "#log2(@-steps), int/double: log2 of number of equilibriation steps, i.e. 6 gives 2^6 steps" << endl;
        cout << "alpha, double: WF parameter for simple gaussian. Analytical sol alpha = omega/2" << endl;
        cout << "stepLenght, double: How far should I move a particle at each MC cycle?" << endl;
        cout << "Importantce sampling?, bool: If the Metropolis Hasting algorithm is used. Then stepLength serves as Delta t" << endl;
        cout << "analytical?, bool: If the analytical expression should be used. Defaults to true" << endl;
        cout << "learningRate *(ln(#particles) + 1)^-1, double: Learning rate for gradient descent. Defaults to 0.01*(ln(#particles) + 1)^-1" << endl;
        cout << "epsilon, double: Tolerance for gradient descent. Defaults to 0.01" << endl;
        cout << "interaction?, bool: If the interacting gaussian should be used. Defaults to false" << endl;
        cout << "filename, string: If the results should be dumped to a file, give the file name. If none is given, a simple print is performed." << endl;
        return 0;
    }

    // Check how many arguments are given and overwrite defaults. Works serially,
    // meaning if 4 parameters are given the first 4 paramters will be
    // overwritten, the rest will be defaults.
    if (argv >= 2)
        numberOfDimensions = (unsigned int)atoi(argc[1]);
    if (argv >= 3)
        numberOfParticles = (unsigned int)atoi(argc[2]);
    if (argv >= 4)
        numberOfMetropolisSteps = (unsigned int)pow(2, atof(argc[3]));
    if (argv >= 5)
        numberOfEquilibrationSteps = (unsigned int)pow(2, atof(argc[4]));
    if (argv >= 6)
        alpha = (double)atof(argc[5]);
    if (argv >= 7)
        stepLength = (double)atof(argc[6]);
    if (argv >= 8)
        importanceSampling = (bool)atoi(argc[7]);
    if (argv >= 9)
        analytical = (bool)atoi(argc[8]);
    if (argv >= 10)
        lr = (double)atof(argc[9]) / (log(numberOfParticles) + 1);
    if (argv >= 11)
        epsilon = (double)atof(argc[10]);
    if (argv >= 12)
        interacting = (bool)atoi(argc[11]);
    if (argv >= 13)
        filename = argc[12];

    // The random engine can also be built without a seed
    auto rng = std::make_unique<Random>(seed);

    // Initialize particles
    auto particles = setupRandomUniformInitialState(
        stepLength, omega, numberOfDimensions, numberOfParticles, *rng);

    // Construct a unique pointer to a new System
    auto hamiltonian = std::make_unique<HarmonicOscillator>(omega);

    // Initialise SimpleGaussian by default
    std::unique_ptr<class WaveFunction> wavefunction;

    if (!analytical)
        wavefunction = std::make_unique<SimpleGaussianNumerical>(alpha, dx); // constructor (can only be moved once).

    if (!interacting)
        wavefunction = std::make_unique<SimpleGaussian>(
            alpha);
    else
    {
        wavefunction = std::make_unique<InteractingGaussian>(
            alpha,
            interactionTerm,
            numberOfParticles);
    }
    // Empty solver pointer, since it uses "rng" in its constructor (can only be
    // moved once).
    std::unique_ptr<class MonteCarlo> solver;

    // Set what solver to use, pass on rng and additional parameters
    if (importanceSampling)
    {
        solver = std::make_unique<MetropolisHastings>(std::move(rng), stepLength, D);
    }
    else
    {
        solver = std::make_unique<Metropolis>(std::move(rng));
    }

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

    // Run the Metropolis algorithm with gradient descent. Notice that the equilibration steps for GD is done inside Optimizer
    auto sampler = system->optimizeMetropolis(
        *system, filename, stepLength, numberOfMetropolisSteps, numberOfEquilibrationSteps, epsilon, lr);

    // Output general information from the simulation, either as file or print
    sampler->output(*system, filename, omega, analytical, importanceSampling);

    return 0;
}