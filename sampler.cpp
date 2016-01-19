#include <iostream>
#include <cmath>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;


Sampler::Sampler(System* system) {
    this->system = system;
    this->stepNumber = 0;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    this->numberOfMetropolisSteps = steps;
}

void Sampler::sample(bool acceptedStep) {
    if (this->stepNumber == 0) {
        this->cumulativeEnergy          = 0;
        this->cumulativeEnergy2         = 0;
        this->cumulativeAcceptanceRate  = 0;
    }

    double localEnergy = this->system->getHamiltonian()->
                         computeLocalEnergy(this->system->getParticles());
    this->cumulativeEnergy  += localEnergy;
    this->cumulativeEnergy2 += localEnergy * localEnergy;
    this->cumulativeAcceptanceRate += acceptedStep;
    this->stepNumber++;
}

void Sampler::printOutputToTerminal() {
    int     np = this->system->getNumberOfParticles();
    int     nd = this->system->getNumberOfDimensions();
    int     ms = this->system->getNumberOfMetropolisSteps();
    int     p  = this->system->getWaveFunction()->getNumberOfParameters();
    double  ef = this->system->getEquilibrationFraction();
    double* pa = this->system->getWaveFunction()->getParameters();

    cout << endl;
    cout << "  -- Physical system info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : " << ms << endl;
    cout << " Number of equilibration steps  : " << std::round(ms*ef) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa[i] << endl;
    }
    cout << endl;
    cout << "  -- Reults -- " << endl;
    cout << " Energy   : " << this->energy << endl;
    cout << " Variance : " << this->energyVariance << endl;
    cout << " Acceptance rate : " << this->acceptanceRate << endl;
    cout << endl;
}

void Sampler::computeAverages() {
    double  equilibrationFraction = this->system->getEquilibrationFraction();
    int     numberOfStepsSampled = this->numberOfMetropolisSteps *
                                   (1.0 - equilibrationFraction);
    double energy2 = this->cumulativeEnergy2 / ((double) numberOfStepsSampled);
    this->energy   = this->cumulativeEnergy  / ((double) numberOfStepsSampled);
    this->energyVariance = energy2 - this->energy * this->energy;
    this->acceptanceRate = this->cumulativeAcceptanceRate / ((double) numberOfStepsSampled);
}
