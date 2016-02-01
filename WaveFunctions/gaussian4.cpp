#include "gaussian4.h"
#include <cmath>
#include "../particle.h"

Gaussian4::Gaussian4(System* system, double alpha) :
        WaveFunction(system) {
    m_numberOfParameters = 1;
    m_parameters = new double[m_numberOfParameters];
    m_parameters[0] = alpha;
}

double Gaussian4::evaluate(Particle* particles) {
    double x = particles[0].getPosition()[0];
    double a = m_parameters[0];
    return std::exp(- a * x*x*x*x);
}

double Gaussian4::computeKineticEnergy(Particle* particles) {
    double x = particles[0].getPosition()[0];
    double a = m_parameters[0];
    return (4 * a * x*x * std::exp(- a * x*x*x*x) * (4 * a * x*x*x*x - 3)) /
            evaluate(particles);
}


