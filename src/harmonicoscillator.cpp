#include <memory>
#include <cassert>
#include <iostream>

#include "../include/harmonicoscillator.h"
#include "../include/particle.h"
#include "../include/wavefunction.h"

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
    auto kineticEnergy = -0.5 * waveFunction.computeLocalLaplasian(particles);
    auto potentialEnergy = 0.0;

    for (auto &particle : particles)
    {
        auto position = particle->getPosition();
        for (size_t pos_index = 0; pos_index < particle->getNumberOfDimensions(); ++pos_index)
            potentialEnergy += position[pos_index] * position[pos_index];
    }
    potentialEnergy *= 0.5 * m_omega * m_omega;

    return kineticEnergy + potentialEnergy;
}
