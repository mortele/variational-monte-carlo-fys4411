#include "multiparticleho.h"
#include <cmath>
#include "../system.h"
#include "../particle.h"


MultiparticleHO::MultiparticleHO(System* system, double alpha) :
        WaveFunction(system) {

    m_numberOfDimensions = system->getNumberOfDimensions();
    m_numberOfParticles = system->getNumberOfParticles();
    m_parameters = new double[1];
    m_parameters[0] = alpha;
}

double MultiparticleHO::evaluate(Particle* particles) {
    double exponentSum = 0;
    for (int i=0; i < m_numberOfParticles; i++) {
        for (int j=0; j < m_numberOfDimensions; j++) {
            exponentSum -= m_parameters[0] * particles[i].getPosition()[j] *
                                             particles[i].getPosition()[j];
        }
    }
    return std::exp(exponentSum);
}

double MultiparticleHO::computeKineticEnergy(Particle* particles) {
    double kineticEnergy = 0;
    for (int i=0; i < m_numberOfParticles; i++) {
        double r2 = 0;
        for (int j=0; j < m_numberOfDimensions; j++) {
            r2 += particles[i].getPosition()[j] *
                  particles[i].getPosition()[j];
        }
        kineticEnergy += m_numberOfDimensions*m_parameters[0] -
                         2*m_parameters[0]*m_parameters[0] * r2;
    }
    return kineticEnergy;
}
