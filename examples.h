#pragma once
#include <iostream>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/gaussian4.h"
#include "WaveFunctions/multiparticleho.h"
#include "WaveFunctions/multiparticlehointeracting.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/harmonicoscillatorinteracting.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"

System* NonInteractingHO(int numberOfDimensions, int numberOfParticles) {
    double omega            = 2.0;
    double alpha            = 0.5*omega;
    double stepLength       = 0.5;
    double equilibrationFraction = 0.0;

    System* system = new System();
    system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setHamiltonian(new HarmonicOscillator(system, omega));
    system->setWaveFunction(new MultiparticleHO(system, alpha));
    system->setEquilibrationFraction(equilibrationFraction);
    system->setStepLength(stepLength);
    system->runMetropolisSteps((int) 1e5);
    return system;
}


System* InteractingHO(int numberOfDimensions, int numberOfParticles) {
    double omega            = 1.0;
    double gamma            = 1.0;
    double alpha            = 0.5*omega;
    double stepLength       = 0.5;
    double equilibrationFraction = 0.0;

    System* system = new System();
    system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setHamiltonian(new HarmonicOscillatorInteracting(system, gamma));
    system->setWaveFunction(new MultiparticleHOInteracting(system, alpha));
    system->setEquilibrationFraction(equilibrationFraction);
    system->setStepLength(stepLength);
    system->runMetropolisSteps((int) 1e2);
    return system;
}
