#include "simplegaussian.h"
#include <cmath>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

SimpleGaussian::SimpleGaussian(System* system, double alpha) :
        WaveFunction(system) {
    this->numberOfParameters = 1;
    this->parameters = new double[this->numberOfParameters];
    this->parameters[0] = alpha;
}

double SimpleGaussian::evaluate(Particle* particle) {
    double position = particle->getPosition()[0];
    return std::exp(-this->parameters[0] * position * position / 2.0);
}

double SimpleGaussian::computeDoubleDerivative(Particle* particles) {
    double position = particles->getPosition()[0];
    return this->parameters[0]*std::exp(-this->parameters[0] * position * position / 2.0) *
            (this->parameters[0] * position * position - 1);
}
