#include "simplegaussian.h"
#include <cmath>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

SimpleGaussian::SimpleGaussian(System* system, double alpha) :
        WaveFunction(system) {
    m_numberOfParameters = 1;
    m_parameters = new double[m_numberOfParameters];
    m_parameters[0] = alpha;
}

double SimpleGaussian::evaluate(Particle* particle) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     */
    return 0;
}

double SimpleGaussian::computeDoubleDerivative(Particle* particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */
    return 0;
}
