#include "system.h"
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"


System::System() {
}

bool System::metropolisStep() {
    int     particle        = Random::nextInt(this->numberOfParticles);
    int     dimension       = Random::nextInt(this->numberOfDimensions);
    double  proposedChange  = (Random::nextDouble()*2-1) * this->stepLength;
    double  waveFunctionOld = this->waveFunction->evaluate(particles);
    this->particles[particle].adjustPosition(proposedChange, dimension);
    double waveFunctionNew  = this->waveFunction->evaluate(particles);

    double waveFunctionSquaredRatio =  waveFunctionNew * waveFunctionNew /
                                      (waveFunctionOld * waveFunctionOld);
    if (waveFunctionSquaredRatio < 1.0) {
        if (waveFunctionSquaredRatio < Random::nextDouble()) {
            this->particles[particle].adjustPosition(-proposedChange, dimension);
            return false;
        }
    } else {
        return true;
    }
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    this->particles                 = this->initialState->getParticles();
    this->sampler                   = new Sampler(this);
    this->numberOfMetropolisSteps   = numberOfMetropolisSteps;
    this->sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    for (int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = this->metropolisStep();

        if (i > this->equilibrationFraction * numberOfMetropolisSteps) {
            sampler->sample(acceptedStep);
        }
    }
    this->sampler->computeAverages();
    this->sampler->printOutputToTerminal();
}

void System::setNumberOfParticles(int numberOfParticles) {
    this->numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    this->numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    this->stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    this->equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    this->hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    this->waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    this->initialState = initialState;
}


