#include <memory>
#include <cmath>
#include <cassert>

#include "simplegaussian.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>

SimpleGaussian::SimpleGaussian(double alpha)
{
    assert(alpha >= 0);
    m_numberOfParameters = 1; // this should not be hard coded
    m_parameters.reserve(1);
    m_parameters.push_back(alpha); // m_parameters is the vector of variational parameters
}

double SimpleGaussian::evaluate(std::vector<std::unique_ptr<class Particle>> &particles)
{
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i]->getPosition()
     * function.
     */
    int num_particles = particles.size();
    int numberOfDimensions = particles.at(0)->getNumberOfDimensions();

    double r2 = 0;  // r2 is the sum of the squared coordinates of the r vector
    double r_q = 0; // r_q is the q'th coordinate of the r vector
    double alpha = m_parameters.at(0);

    for (int i = 0; i < num_particles; i++)
    {
        Particle particle = *particles.at(i);
        for (int q = 0; q < numberOfDimensions; q++)
        {
            r_q = particle.getPosition().at(q);
            r2 += r_q * r_q;
        }
    }

    return std::exp(-alpha * r2);
}

double SimpleGaussian::computeParamDerivative(std::vector<std::unique_ptr<class Particle>> &particles, int parameterIndex)
{
    /* Note that by derivative, we actually
     * mean the derivative with respect to the variational parameters.
     */
    int num_particles = particles.size();
    int numberOfDimensions = particles.at(0)->getNumberOfDimensions();
    double parameter = m_parameters.at(parameterIndex); // this is not used now, but can be used when we generalize

    double r2_sum = 0;
    double r_q = 0;

    for (int k = 0; k < num_particles; k++)
    {
        Particle particle = *particles.at(k);
        for (int q = 0; q < numberOfDimensions; q++)
        {
            r_q = particle.getPosition().at(q);
            r2_sum += r_q * r_q;
        }
    }
    // notice this does not depend on the param as it is divided by the WF.
    return -r2_sum; // analytic derivative wrt alpha, only 1 param to optimize now, This needs to be generalized
}

double SimpleGaussian::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>> &particles)
{
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */
    int num_particles = particles.size();
    int numberOfDimensions = particles.at(0)->getNumberOfDimensions();
    double alpha = m_parameters.at(0);

    double r2_sum = 0;
    double r_q = 0;

    for (int k = 0; k < num_particles; k++)
    {
        Particle particle = *particles.at(k);
        for (int q = 0; q < numberOfDimensions; q++)
        {
            r_q = particle.getPosition().at(q);
            r2_sum += r_q * r_q;
        }
    }

    return 2 * alpha * (2 * alpha * r2_sum - num_particles * numberOfDimensions); // analytic double derivative
}

double SimpleGaussian::evaluate_w(int proposed_particle_idx, class Particle &proposed_particle, class Particle &old_particle, std::vector<std::unique_ptr<class Particle>> &particles)
{
    /*
     This is the wave function ratio for the Metropolis algorithm.
     It is a clever way to avoid having to evaluate the wave function for all particles at each step.
     */
    static const int numberOfDimensions = particles.at(0)->getNumberOfDimensions(); // static to avoid redeclaration between calls
    double alpha = m_parameters.at(0);

    double r2_proposed, r2_old;
    r2_proposed = 0;
    r2_old = 0;

    for (int i = 0; i < numberOfDimensions; i++)
    {
        r2_proposed += proposed_particle.getPosition().at(i) * proposed_particle.getPosition().at(i);
        r2_old += old_particle.getPosition().at(i) * old_particle.getPosition().at(i);
    }

    return std::exp(-2.0 * alpha * (r2_proposed - r2_old));
}

void SimpleGaussian::quantumForce(std::vector<std::unique_ptr<class Particle>> &particles, Particle &particle, std::vector<double> &force)
{
    static const int numberOfDimensions = particle.getNumberOfDimensions(); // static to avoid redeclaration between calls
    double alpha = m_parameters.at(0);

    for (int q = 0; q < numberOfDimensions; q++)
    {
        force.at(q) = -4.0 * alpha * particle.getPosition().at(q);
    }
}

SimpleGaussianNumerical::SimpleGaussianNumerical(double alpha, double dx) : SimpleGaussian(alpha)
{
    m_dx = dx;
}

double SimpleGaussianNumerical::evaluateSingleParticle(class Particle &particle) 
{
    double r2 = 0;
    static const int numberOfDimensions = particle.getNumberOfDimensions();

    for(int q = 0; q < numberOfDimensions; q++) {
        r2 += particle.getPosition().at(q)*particle.getPosition().at(q);
    }

    return std::exp(-m_parameters.at(0)*r2);
}

double SimpleGaussianNumerical::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>> &particles)
{
    int num_particles = particles.size();
    int numberOfDimensions = particles.at(0)->getNumberOfDimensions();
    double der_sum = 0;
    double r_q, gx, gxpdx, gxmdx, der;

    for (int i = 0; i < num_particles; i++)
    {
        Particle &particle = *particles.at(i);

        gx = evaluateSingleParticle(particle);          // gx = g(x) this is the value of the wave function at the current position
        for (int q = 0; q < numberOfDimensions; q++)
        {
            r_q = particle.getPosition().at(q);
            particle.adjustPosition(m_dx, q);               // adjust the position of the particle by dx in the qth dimension
            gxpdx = evaluateSingleParticle(particle);                    // gxpdx = g(x + dx) this is the value of the wave function at the new position
            particle.adjustPosition(-2 * m_dx, q);          // adjust the position of the particle by -2dx in the qth dimension
            gxmdx = evaluateSingleParticle(particle);                    // gxmdx = g(x - dx) this is the value of the wave function at the new position
            der = (gxpdx - 2 * gx + gxmdx) / (m_dx * m_dx); // der = double derivative
            particle.setPosition(r_q, q);                   // reset the position of the particle to the original value
            der_sum += der/gx;                                 // sum the double derivatives over all particles and dimensions
        }
    }

    return der_sum; // divide by the value of the wave function at the current position
}

void SimpleGaussian::setParameters(std::vector<double> parameters)
{
    assert((int)parameters.size() == m_numberOfParameters);
    m_parameters = parameters;
}