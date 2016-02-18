#pragma once
#include "wavefunction.h"

class MultiparticleHO : public WaveFunction {
public:
    MultiparticleHO(class System* system, double alpha);
    virtual double evaluate(class Particle* particles);
    virtual double computeKineticEnergy(class Particle* particles);

private:
    int m_numberOfDimensions = 0;
    int m_numberOfParticles = 0;
};

