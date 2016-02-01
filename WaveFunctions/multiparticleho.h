#pragma once
#include "wavefunction.h"

class MultiparticleHO : public WaveFunction {
public:
    MultiparticleHO(class System* system, double alpha);
    double evaluate(class Particle* particles);
    double computeDoubleDerivative(class Particle* particles);

private:
    int m_numberOfDimensions = 0;
    int m_numberOfParticles = 0;
};

