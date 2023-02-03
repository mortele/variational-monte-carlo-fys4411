#include <memory>
#include <cmath>
#include <cassert>

#include "simplegaussian.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include<iostream>

SimpleGaussian::SimpleGaussian(double alpha)
{
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(std::vector<std::unique_ptr<class Particle>>& particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i]->getPosition()
     * function.
     */

    double psi_T = 1;
    int num_particles = particles.size();
    int numberOfDimensions = particles.at(0)->getNumberOfDimensions();
    
    double r2 = 0;
    double r_q = 0;
    double alpha = m_parameters.at(0);

    for(int i = 0; i < num_particles; i++) {
        Particle& particle = *particles.at(i);
        r2 = 0;
        for(int q = 0; q < numberOfDimensions; q++) {
            r_q =  particle.getPosition().at(q);
            r2 += r_q * r_q;
        }

        psi_T *= std::exp( -alpha * r2 );
    }

    return psi_T;
}

double SimpleGaussian::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles) {
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
    
    // double psi_T = evaluate(particles);
    double sum = 0;

    double r2_sum = 0;
    double r_q = 0;

    for(int k = 0; k < num_particles; k++) {
        Particle& particle = *particles.at(k);
        for(int q = 0; q < numberOfDimensions; q++) {
            r_q = particle.getPosition().at(q);
            r2_sum += r_q * r_q;
        }
    }

    return 2*alpha*(2*alpha*r2_sum - num_particles*numberOfDimensions);
}