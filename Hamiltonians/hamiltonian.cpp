#include "hamiltonian.h"
#include <iostream>
#include "system.h"
#include "WaveFunctions/wavefunction.h"
#include "particle.h"

using std::cout;
using std::endl;


Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

void Hamiltonian::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

double Hamiltonian::getExactEnergy() {
    return m_exactEnergy;
}

double Hamiltonian::computeKineticEnergy(Particle* particles) {
    const double WaveFunction = 2*m_waveFunction->evaluate(particles);
    const double dx = 1e-5;
    double doubleDerivative = 0;

    for (int i=0; i<m_system->getNumberOfParticles(); i++) {
        for (int j=0; j<m_system->getNumberOfDimensions(); j++) {
            particles[i].adjustPosition(dx, j);
            const double waveFunctionPlus = m_waveFunction->evaluate(particles);
            particles[i].adjustPosition(-2*dx, j);
            const double waveFunctionMinus = m_waveFunction->evaluate(particles);
            particles[i].adjustPosition(dx, j);
            doubleDerivative += waveFunctionPlus - WaveFunction + waveFunctionMinus;
        }
    }
    doubleDerivative /= (dx*dx);
    return -doubleDerivative / WaveFunction;
}

