#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

SimpleGaussian::SimpleGaussian(System* system, double alpha, double dt) :
        WaveFunction(system) {
    assert(alpha >= 0);
    assert(dt >= 0);
    m_numberOfParameters = 2;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(dt);
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles, int particle) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */
	double r2 = 0;
	for( int dim = 0; dim < m_system->getNumberOfDimensions(); dim++ )
	{
		r2 += std::pow(particles[particle]->getPosition()[dim], 2);
	}
	//return -m_parameters[0]*r2;
	return std::exp(-m_parameters[0]*r2);
}

double SimpleGaussian::computeDerivative(std::vector<class Particle*> particles) {
    /* Computes the derivative of wave function ansatz as function of variational parameters.
     * (I think this is correct)
     */
	double r2 = 0;
    for (int particle = 0; particle < m_system->getNumberOfParticles(); particle++)
    {
        for( int dim = 0; dim < m_system->getNumberOfDimensions(); dim++ )
        {
            r2 += std::pow(particles[particle]->getPosition()[dim], 2);
        }
    }
	return -.5*r2;
}


double SimpleGaussian::computeDoubleDerivative(double r2) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */

	/* double r2 = 0;
	for( int dim = 0; dim < m_system->getNumberOfDimensions(); dim++ )
	{
		r2 += std::pow(particles[particle]->getPosition()[dim], 2);
	} */

    // since HarmonicOscillator::computeLocalEnergy allready finds r2, we can just reuse
    // that one. This might have to change in the future though, if we need the double
    // derivative somewhere else
	double D = m_system->getNumberOfDimensions();
	double N = m_system->getNumberOfParticles();
	return -2*D*N*m_parameters[0] + 4*pow(m_parameters[0], 2)*r2;
}





































