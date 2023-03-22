#pragma once

#include <vector>
#include <memory>

/// @brief Base class for Monte Carlo solvers
class MonteCarlo
{
public:
    /// @brief Constructor for MonteCarlo
    /// @param rng random number generator
    MonteCarlo(std::unique_ptr<class Random> rng);
    /// @brief Destructor for MonteCarlo
    virtual ~MonteCarlo() = default;

    /// @brief Perform a Monte Carlo step
    /// @param stepLength step length
    /// @param waveFunction wave function
    /// @param particles vector of particles
    /// @return true if step was accepted
    virtual bool step(
        double stepLength,
        class WaveFunction &waveFunction,
        std::vector<std::unique_ptr<class Particle>> &particles) = 0;

protected:
    std::unique_ptr<class Random> m_rng;
};
