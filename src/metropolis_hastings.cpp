#include <memory>
#include <vector>

#include "../include/metropolis_hastings.h"
#include "../include/wavefunction.h"
#include "../include/particle.h"
#include "../include/random.h"

MetropolisHastings::MetropolisHastings(std::unique_ptr<class Random> rng)
    : MonteCarlo(std::move(rng))
{
}

bool MetropolisHastings::step(
    double stepLength,
    class WaveFunction &waveFunction,
    std::vector<std::unique_ptr<class Particle>> &particles)
{
    // Choose a random particle
    auto particle_index = m_rng->nextInt(particles.size() - 1);
    auto particle = *particles[particle_index];

    // Save old particle
    auto new_particle = particle;

    // Propose a new position
    for (size_t pos_index = 0; pos_index < particle.getNumberOfDimensions(); ++pos_index)
        new_particle.adjustPosition(m_rng->nextGaussian(0, stepLength), pos_index);

    // copy particles array 
    auto new_particles = std::vector<std::unique_ptr<class Particle>>(particles.size());
    for (size_t i = 0; i < particles.size(); ++i)
        new_particles[i] = std::make_unique<Particle>(*particles[i]);

    new_particles[particle_index] = std::make_unique<Particle>(new_particle);

    // Calculate the ratio of the new and old wave functions
    auto ratio = waveFunction.evaluateRatio(new_particles, particles);

    // Accept or reject the move
    if (m_rng->nextDouble() < ratio * ratio)
    {
        // Accept the move
        particles[particle_index] = std::make_unique<Particle>(new_particle);
        return true;
    }

    return false;
}
