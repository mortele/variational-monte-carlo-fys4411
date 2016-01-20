#include "hamiltonian.h"
#include "../system.h"

Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
    m_waveFunction = system->getWaveFunction();
}

void Hamiltonian::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

