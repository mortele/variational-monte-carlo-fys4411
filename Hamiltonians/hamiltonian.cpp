#include "hamiltonian.h"

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

}

