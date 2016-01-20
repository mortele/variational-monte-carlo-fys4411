#include "randomuniform.h"
#include "Math/random.h"
#include "../particle.h"
#include "../system.h"


RandomUniform::RandomUniform(System*    system,
                             int        numberOfDimensions,
                             int        numberOfParticles)  :
        InitialState(system) {
    m_numberOfDimensions = numberOfDimensions;
    m_numberOfParticles  = numberOfParticles;
    m_system->setNumberOfDimensions(numberOfDimensions);
    m_system->setNumberOfParticles(numberOfParticles);
    setupInitialState();
}

void RandomUniform::setupInitialState() {
    m_particles = new Particle[m_numberOfParticles];

    for (int i=0; i < m_numberOfParticles; i++) {
        double* position = new double[m_numberOfDimensions];
        for (int j=0; j < m_numberOfDimensions; j++) {
            position[j] = Random::nextDouble()*2-1;
        }
        m_particles[i].setNumberOfDimensions(m_numberOfDimensions);
        m_particles[i].setPosition(position);
    }
}
