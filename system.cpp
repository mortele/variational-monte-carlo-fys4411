#include "system.h"
#include <iostream>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"

using std::cout;
using std::endl;

bool System::metropolisStep() {
    int     particle        = Random::nextInt(m_numberOfParticles);
    int     dimension       = Random::nextInt(m_numberOfDimensions);
    double  proposedChange  = (Random::nextDouble()*2-1) * m_stepLength;

    double  waveFunctionOld = m_waveFunction->evaluate(m_particles);
    m_particles[particle].adjustPosition(proposedChange, dimension);
    double waveFunctionNew  = m_waveFunction->evaluate(m_particles);

    double waveFunctionSquaredRatio =  waveFunctionNew * waveFunctionNew /
                                      (waveFunctionOld * waveFunctionOld);
    if (waveFunctionSquaredRatio < 1.0) {
        if (waveFunctionSquaredRatio < Random::nextDouble()) {
            m_particles[particle].adjustPosition(-proposedChange, dimension);
            return false;
        }
    }
    return true;
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    m_hamiltonian->setWaveFunction(m_waveFunction);

    for (int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = metropolisStep();

        if (i > m_equilibrationFraction * numberOfMetropolisSteps) {
            m_sampler->sample(acceptedStep);
        }
    }
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}


