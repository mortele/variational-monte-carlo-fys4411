#include <memory>
#include <vector>

#include "metropolis.h"
#include "WaveFunctions/wavefunction.h"
#include "particle.h"
#include "Math/random.h"

#include <iostream>

Metropolis::Metropolis(std::unique_ptr<class Random> rng)
    : MonteCarlo(std::move(rng))
{
}

bool Metropolis::step(
    double stepLength,
    class WaveFunction &waveFunction,
    std::vector<std::unique_ptr<class Particle>> &particles)
{
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change its position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated at
     * this new position with the one at the old position).
     */

    int numberOfParticles = particles.size();
    int numberOfDimensions = particles.at(0)->getNumberOfDimensions();

    // double Psi_old = waveFunction.evaluate(particles);

    int proposed_particle_idx = m_rng->nextInt(0, numberOfParticles - 1); // Choose a particle at random
    Particle &proposed_particle = *particles.at(proposed_particle_idx);   // Get a reference to the particle
    Particle old_particle = proposed_particle;                            // Save the old position

    for (int q = 0; q < numberOfDimensions; q++)
        proposed_particle.adjustPosition(stepLength * (m_rng->nextDouble() - .5), q);

    // double Psi_new = waveFunction.evaluate(particles);

    // double w = (Psi_new * Psi_new) / (Psi_old * Psi_old); // Metropolis test
    double w = waveFunction.evaluate_w(proposed_particle_idx,
                                       proposed_particle,
                                       old_particle,
                                       particles);

    if (w >= m_rng->nextDouble()) // Accept the step
    {
        return true;
    }
    else
    {
        *particles.at(proposed_particle_idx) = old_particle; // Reject the step and restore the old position
        return false;
    }
}