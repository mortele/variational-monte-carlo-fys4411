#include "heliumatom.h"
#include <cmath>
#include "particle.h"
#include "system.h"

HeliumAtom::HeliumAtom(System* system, bool interaction) :
        Hamiltonian(system) {
    m_interaction = interaction;
}

HeliumAtom::HeliumAtom(System* system) :
        HeliumAtom(system, true) {
}

double HeliumAtom::computeLocalEnergy(Particle* particles) {
    double r12  = 0;
    double r1   = 0;
    double r2   = 0;
    for (int k=0; k<m_system->getNumberOfDimensions(); k++) {
        const double x1     = particles[0].getPosition()[k];
        const double x2     = particles[1].getPosition()[k];
        const double x12    = x2 - x1;
        r1  += x1  * x1;
        r2  += x2  * x2;
        r12 += x12 * x12;
    }
    const double kineticEnergy      = Hamiltonian::computeKineticEnergy(particles);
    const double interactionEnergy  = (m_interaction ? 1.0 / sqrt(r12) : 0.0);
    const double potentialEnergy    = -2.0 / sqrt(r1) - 2.0 / sqrt(r2);
    return kineticEnergy + potentialEnergy + interactionEnergy;
}

