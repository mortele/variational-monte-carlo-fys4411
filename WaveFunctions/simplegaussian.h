#pragma once
#include <memory>
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(std::shared_ptr<class System> system, double alpha);
    double evaluate(std::vector<std::unique_ptr<class Particle>> particles);
    double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>> particles);
};
