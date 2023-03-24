#pragma once
#include <memory>
#include <vector>

#include "hamiltonian.h"
#include "particle.h"

class HarmonicOscillator : public Hamiltonian
{
public:
    /// @brief Constructor for the HarmonicOscillator class.
    HarmonicOscillator(double omega);
    /// @brief Compute the local energy of the system.
    /// @param waveFunction The trial wave function to use.
    /// @param particles Vector of particles.
    /// @return The local energy in natural units with m = 1.
    double computeLocalEnergy(
        class WaveFunction &waveFunction,
        std::vector<std::unique_ptr<class Particle>> &particles);

private:
    /// @brief Omega frequency parameter in harmonic oscillator.
    double m_omega;
};
