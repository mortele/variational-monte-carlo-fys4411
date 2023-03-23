#include "wavefunction.h"
#include "../particle.h"

double WaveFunction::computeLocalLaplasian(std::vector<std::unique_ptr<class Particle>> &particles)
{
    // Compute the local energy by numerical differentiation
    double h = 1e-4; // never do this, but I assume dimensionless units, so fine
    double laplasian = 0;

    for (auto &particle : particles)
    {
        for (size_t i = 0; i < particle->getPosition().size(); i++)
        {
            auto mid_val = evaluate(particles);
            particle->adjustPosition(h, i);
            auto plus_val = evaluate(particles);
            particle->adjustPosition(-2 * h, i);
            auto minus_val = evaluate(particles);
            particle->adjustPosition(h, i);
            laplasian += (plus_val - 2 * mid_val + minus_val) / (h * h);
        }
    }
    return laplasian / evaluate(particles);
}

double WaveFunction::evaluateRatio(std::vector<std::unique_ptr<class Particle>> &particles_numerator, std::vector<std::unique_ptr<class Particle>> &particles_denominator)
{
    return evaluate(particles_numerator) / evaluate(particles_denominator);
}