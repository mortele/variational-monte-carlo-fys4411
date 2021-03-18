#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha, double dt);
    double evaluate(std::vector<class Particle*> particles);
    double computeDerivative(std::vector<class Particle*> particles);
    double computeDoubleDerivative(double r2);
    // double computeDoubleDerivative(std::vector<class Particle*> particles, int particle);
};
