#include <memory>
#include <cmath>
#include <cassert>

#include "../include/simplegaussian.h"
#include "../include/particle.h"

SimpleGaussian::SimpleGaussian(double alpha)
{
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(std::vector<std::unique_ptr<class Particle>> &particles)
{
    double psi = 1.0;
    double alpha = m_parameters[0]; // alpha is the first and only parameter for now.

    for (size_t i = 0; i < particles.size(); i++)
    {
        // Let's support as many dimensions as we want.
        double r2 = 0;
        for (size_t j = 0; j < particles[i]->getPosition().size(); j++)
            r2 += particles[i]->getPosition()[j] * particles[i]->getPosition()[j];
        // spherical ansatz
        double g = exp(-alpha * r2);

        // Trial wave function is product of g for all particles.
        // f ignored for now, due to considering non interacting particles.
        psi = psi * g;
    }
    return psi;
}


double SimpleGaussian::computeLocalLaplasian(std::vector<std::unique_ptr<class Particle>> &particles)
{
    // The expression I got for a single laplasian is, in invariant form, follows:
    // (4 * alpha^2 * r_i^2 - 2 * alpha * NDIM)
    // so it takes to sum over all particles.
    double alpha = m_parameters[0];
    double sum_laplasian = 0.0;
    for (size_t i = 0; i < particles.size(); i++)
    {
        double r2 = 0.0;
        for (size_t j = 0; j < particles[i]->getPosition().size(); ++j)
            r2 += particles[i]->getPosition()[j] * particles[i]->getPosition()[j];
        sum_laplasian += 4 * alpha * alpha * r2 - 2 * alpha * particles[i]->getPosition().size();
    }
    return sum_laplasian;
}

double SimpleGaussian::evaluateRatio(std::vector<std::unique_ptr<class Particle>> &particles_numerator, std::vector<std::unique_ptr<class Particle>> &particles_denominator)
{
    assert(particles_numerator.size() == particles_denominator.size());
    double ratio = 1.0;
    double alpha = m_parameters[0];

    for (size_t i = 0; i < particles_numerator.size(); i++)
    {
        double r2_numerator = 0.0;
        double r2_denominator = 0.0;
        for (size_t j = 0; j < particles_numerator[i]->getPosition().size(); j++)
        {
            r2_numerator += particles_numerator[i]->getPosition()[j] * particles_numerator[i]->getPosition()[j];
            r2_denominator += particles_denominator[i]->getPosition()[j] * particles_denominator[i]->getPosition()[j];
        }
        ratio *= exp(-alpha * (r2_numerator - r2_denominator));
    }
    return ratio;
}
