#include "harmonicoscillator.h"
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
        Hamiltonian(system) {
    this->omega  = omega;
    this->omega2 = omega*omega;
}

double HarmonicOscillator::computeLocalEnergy(Particle* particles) {
    double potentialEnergy = 0;
    for (int i=0; i < this->system->getNumberOfParticles(); i++) {
        double* position = particles[i].getPosition();
        double  r2       = 0;

        for (int j=0; j < this->system->getNumberOfDimensions(); j++) {
            r2 += position[j]*position[j];
        }
        potentialEnergy += 0.5 * this->omega2 * r2;
    }
    WaveFunction* waveFunction = this->system->getWaveFunction();
    double kineticEnergy = - waveFunction->computeDoubleDerivative(particles) /
                             (2 * waveFunction->evaluate(particles));
    return kineticEnergy + potentialEnergy;
}

