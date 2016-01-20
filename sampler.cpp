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
    m_system = system;
    m_stepNumber = 0;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

void Sampler::sample(bool acceptedStep) {
    if (m_stepNumber == 0) {
        m_cumulativeEnergy          = 0;
        m_cumulativeEnergy2         = 0;
        m_cumulativeAcceptanceRate  = 0;
    }

    double localEnergy = m_system->getHamiltonian()->
                         computeLocalEnergy(m_system->getParticles());
    m_cumulativeEnergy  += localEnergy;
    m_cumulativeEnergy2 += localEnergy * localEnergy;
    m_cumulativeAcceptanceRate += acceptedStep;
    m_stepNumber++;
}

void Sampler::printOutputToTerminal() {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibrationFraction();
    double* pa = m_system->getWaveFunction()->getParameters();
    bool    ek = m_system->getHamiltonian()->getExactGroundStateEnergyKnown();
    double  ee = m_system->getHamiltonian()->getExactEnergy();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa[i] << endl;
    }
    cout << endl;
    cout << "  -- Reults -- " << endl;
    if (ek) {
        cout << " Exact energy    : " << ee << endl;
    }
    cout << " Energy          : " << m_energy << endl;
    cout << " Variance        : " << m_energyVariance << endl;
    cout << " Acceptance rate : " << m_acceptanceRate << endl;
    cout << endl;
}

void Sampler::computeAverages() {
    double  equilibrationFraction = m_system->getEquilibrationFraction();
    int     numberOfStepsSampled = m_numberOfMetropolisSteps *
                                   (1.0 - equilibrationFraction) - 1;
    double energy2 = m_cumulativeEnergy2 / ((double) numberOfStepsSampled);
    m_energy         = m_cumulativeEnergy  / ((double) numberOfStepsSampled);
    m_energyVariance = energy2 - m_energy * m_energy;
    m_acceptanceRate = m_cumulativeAcceptanceRate / ((double) numberOfStepsSampled);
}
