#include <iostream>
#include <fstream>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"

#include <filesystem>
namespace fs = std::filesystem;

using namespace std;


int main() {
    // Seed for the random number generator
    int seed = 2020;

    int numberOfDimensions[]    = {1, 2};
    int numberOfParticles[]     = {1, 2}; //{1,10,100,500};
    int numberOfSteps           = (int) 1e6;
    double omega                = 1.0;              // Oscillator frequency.
    double alpha[]              = {.4};// Variational parameter.
	//double alpha[] = {0.3, 0.34, 0.38, 0.42, 0.46, 0.5, 0.54, 0.58, 0.62, 0.66, 0.7};
    double stepLength           = 2;              // Metropolis step length.
    double equilibration        = 0.1;              // Amount of the total steps used
    // for equilibration.

    // clears output file
    ofstream outfile;
    outfile.open ("res/results.csv", ios::out | ios::trunc);
    outfile << 
    "nParticles;nDimensions;nMetropolisSteps;EquilibrationFraction;foundEnergy;elapsedTime;nParameters;Parameters(undefinedNumber)\n";
    outfile.close();
    outfile.open ("res/energies.csv", ios::out | ios::trunc);
    outfile.close();

    for (int nDim : numberOfDimensions)
    {
        for (int nPar : numberOfParticles)
        {
            for (double nAlpha : alpha)
            {
                System* system = new System(seed);
                system->setHamiltonian              (new HarmonicOscillator(system, omega, false));
                system->setWaveFunction             (new SimpleGaussian(system, nAlpha));
                system->setInitialState             (new RandomUniform(system, nDim, nPar));
                system->setEquilibrationFraction    (equilibration);
                system->setStepLength               (stepLength);
                system->runMetropolisSteps          (numberOfSteps);
            }
        }
    }
    
    return 0;
}
