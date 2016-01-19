#include "hamiltonian.h"

Hamiltonian::Hamiltonian(System* system) {
    this->system = system;
}

void Hamiltonian::setWaveFunction(WaveFunction* waveFunction) {
    this->waveFunction = waveFunction;
}

double Hamiltonian::computeKineticEnergy(Particle* particles) {

}

