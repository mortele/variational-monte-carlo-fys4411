#pragma once

#include <memory>

#include "montecarlo.h"
#include "particle.h"

/// @brief Metropolis Monte Carlo solver
class Metropolis : public MonteCarlo
{
public:
    /// @brief Constructor for Metropolis
    Metropolis(std::unique_ptr<class Random> rng);
    /// @brief Perform a Metropolis Monte Carlo step
    /// @param stepLength length of the step
    /// @param waveFunction wave function to change
    /// @param particles vector of particles
    /// @return true if the step was accepted, false otherwise
    bool step(
        double stepLength,
        class WaveFunction &waveFunction,
        std::vector<std::unique_ptr<class Particle>> &particles);
};
