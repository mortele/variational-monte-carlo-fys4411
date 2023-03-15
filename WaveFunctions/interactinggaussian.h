#pragma once

#include <memory>

#include "wavefunction.h"

class InteractingGaussian : public WaveFunction
{
public:
    InteractingGaussian(double alpha, double a, int num_particles); // a is the interaction parameter
    double evaluate(std::vector<std::unique_ptr<class Particle>> &particles);
    double evaluate_w(int proposed_particle_idx, class Particle &proposed_particle, class Particle &old_particle, std::vector<std::unique_ptr<class Particle>> &particles);
    double computeParamDerivative(std::vector<std::unique_ptr<class Particle>> &particles, int parameterIndex);
    double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>> &particles);
    void quantumForce(Particle &particle, std::vector<double> &force);
    void setParameters(std::vector<double> parameters);
};

class InteractingGaussianNumerical : public InteractingGaussian
{
public:
    InteractingGaussianNumerical(double alpha, double dx, double a, int num_particles);
    double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>> &particles);

private:
    double m_dx;
};