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

using namespace std;


int main() {
    // Seed for the random number generator
    int seed = 2020;

    int numberOfDimensions[]    = {1};
    int numberOfParticles[]     = {1,2,3}; //{1,10,100,500};
    int numberOfSteps           = (int) 1e6;
    double omega                = 1.0;          // Oscillator frequency.
    double alpha                = 0.5;          // Variational parameter.
    double stepLength           = 0.1;          // Metropolis step length.
    double equilibration        = 0.1;          // Amount of the total steps used
    // for equilibration.

    // clears output file
    ofstream outfile;
    outfile.open ("results.txt", ios::out | ios::trunc);
    outfile.close();

    for (int nDim : numberOfDimensions)
    {
            for (int nPar : numberOfParticles)
        {
            System* system = new System(seed);
            system->setHamiltonian              (new HarmonicOscillator(system, omega));
            system->setWaveFunction             (new SimpleGaussian(system, alpha));
            system->setInitialState             (new RandomUniform(system, nDim, nPar));
            system->setEquilibrationFraction    (equilibration);
            system->setStepLength               (stepLength);
            system->runMetropolisSteps          (numberOfSteps);
        }
    }
    
    
    
    
    return 0;
}
