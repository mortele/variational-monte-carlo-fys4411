#include "simplegaussian.h"
#include <cmath>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

SimpleGaussian::SimpleGaussian(System* system, double alpha) :
        WaveFunction(system) {
    m_numberOfParameters = 1;
    m_parameters = new double[m_numberOfParameters];
    m_parameters[0] = alpha;
}

double SimpleGaussian::evaluate(Particle* particle) {
    double x = particle->getPosition()[0];
    double a = m_parameters[0];
    return std::exp(- a * x*x / 2.0);
}

double SimpleGaussian::computeKineticEnergy(Particle* particles) {
    double x = particles->getPosition()[0];
    double a = m_parameters[0];
    return (a * std::exp(- a * x*x / 2.0) * (a * x*x - 1)) / evaluate(particles);
}
