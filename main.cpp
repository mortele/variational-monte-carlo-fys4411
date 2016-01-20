#include <iostream>
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
    System* system = new System();
    system->setHamiltonian(new HarmonicOscillator(system, 1));
    system->setWaveFunction(new SimpleGaussian(system, 1));
    system->setInitialState(new RandomUniform(system, 1, 1));
    system->setEquilibrationFraction(0.1);
    system->setStepLength(0.1);
    system->runMetropolisSteps(1e6);
    return 0;
}

