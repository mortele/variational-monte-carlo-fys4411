#include "harmonicoscillator.h"
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
        Hamiltonian(system) {
    m_omega  = omega;
}

double HarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {
    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * m_waveFunction variable in the super-class in order to compute the
     * kinetic energy here.
     */

    double potentialEnergy = 0;
    double kineticEnergy   = 0;
    return kineticEnergy + potentialEnergy;
}

