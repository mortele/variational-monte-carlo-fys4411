#include "multiparticleho.h"
#include <cmath>
#include "../system.h"
#include "../particle.h"


MultiparticleHO::MultiparticleHO(System* system, double alpha) :
        MultiparticleHO(system, alpha, 1.0) {
    m_numberOfParameters = 1;
}

MultiparticleHO::MultiparticleHO(System* system, double alpha, double beta) :
        WaveFunction(system) {

    m_numberOfParameters = 2;
    m_numberOfDimensions = system->getNumberOfDimensions();
    m_numberOfParticles = system->getNumberOfParticles();
    m_parameters = new double[2];
    m_parameters[0] = alpha;
    m_parameters[1] = beta;
}

double MultiparticleHO::evaluate(Particle* particles) {
    double exponentSum = 0;
    for (int i=0; i < m_numberOfParticles; i++) {
        for (int j=0; j < m_numberOfDimensions; j++) {
            exponentSum +=  particles[i].getPosition()[j] *
                            particles[i].getPosition()[j] *
                            (j==2 ? m_parameters[1] : 1.0);

        }
    }
    return std::exp(- m_parameters[0] * exponentSum);
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
