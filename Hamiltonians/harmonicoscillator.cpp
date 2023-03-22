#include <memory>
#include <cassert>
#include <iostream>

#include "harmonicoscillator.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(double omega)
{
    assert(omega > 0);
    m_omega = omega;
}

double HarmonicOscillator::computeLocalEnergy(
    class WaveFunction &waveFunction,
    std::vector<std::unique_ptr<class Particle>> &particles)
{
    /* Here, you need to compute the kinetic and potential energies.
     * Access to the wave function methods can be done using the dot notation
     * for references, e.g., wavefunction.computeDoubleDerivative(particles),
     * to get the Laplacian of the wave function.
     * */

    double potentialEnergy = 0;
    double kineticEnergy = 0;
    return kineticEnergy + potentialEnergy;
}
