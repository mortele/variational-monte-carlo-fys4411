#pragma once

#include <memory>

#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(double alpha);
    double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
    double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
    void quantumForce(Particle& particle, std::vector<double>& force);
};

class SimpleGaussianNumerical : public SimpleGaussian {
public:
    SimpleGaussianNumerical(double alpha, double dx);
    double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
private:
    double m_dx;

};