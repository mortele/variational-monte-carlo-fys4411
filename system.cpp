#include <memory>
#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"

#include <iostream>


System::System() {
    m_random = std::make_unique<Random>();
}

System::System(int seed) {
    // m_random = std::unique_ptr<Random>(new Random(seed));
    m_random = std::make_unique<Random>(seed);
}

bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */

    return false;
}

void System::runMetropolisSteps(unsigned int numberOfMetropolisSteps) {
    m_particles = m_initialState->getParticles();
    m_sampler = std::make_unique<Sampler>(this);
    m_numberOfMetropolisSteps = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    for (unsigned int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = metropolisStep();

        /* Here you should sample the energy (and maybe other things using
         * the m_sampler instance of the Sampler class. Make sure, though,
         * to only begin sampling after you have let the system equilibrate
         * for a while. You may handle this using the fraction of steps which
         * are equilibration steps; m_equilibrationFraction.
         */
        m_sampler->sample(acceptedStep);
    }
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();
}

void System::setNumberOfParticles(unsigned int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(unsigned int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(std::unique_ptr<class Hamiltonian> hamiltonian) {
    m_hamiltonian = std::move(hamiltonian);
}

void System::setWaveFunction(std::unique_ptr<class WaveFunction> waveFunction) {
    m_waveFunction = std::move(waveFunction);
}

void System::setInitialState(std::unique_ptr<class InitialState> initialState) {
    m_initialState = std::move(initialState);
}


