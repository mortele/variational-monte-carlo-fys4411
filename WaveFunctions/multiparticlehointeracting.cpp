#include "multiparticlehointeracting.h"
#include <cmath>
#include "Hamiltonians/hamiltonian.h"
#include "system.h"
#include "particle.h"
#include <iostream>

using std::cout;
using std::endl;


MultiparticleHOInteracting::MultiparticleHOInteracting(System* system,
                                                       double alpha,
                                                       double beta) :
        MultiparticleHO(system, alpha, beta) {
    m_a     = 0.0043;
    m_a2    = m_a*m_a;
}


double MultiparticleHOInteracting::evaluate(Particle* particles) {
    double nonInteractingWavefunction = MultiparticleHO::evaluate(particles);
    double interactionTerm = 1.0;

    for (int i=0; i<m_system->getNumberOfParticles(); i++) {
        for (int j=i+1; j<m_system->getNumberOfParticles(); j++) {
            double r2 = 0;
            for (int k=0; k<m_system->getNumberOfDimensions(); k++) {
                const double x = particles[i].getPosition()[k] - particles[j].getPosition()[k];
                r2 += x*x;
            }
            const double r = sqrt(r2);
            interactionTerm *= (1.0 - m_a/r) * (r > m_a);
        }
    }
    return nonInteractingWavefunction * interactionTerm;
}
