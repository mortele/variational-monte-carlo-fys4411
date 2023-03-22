#pragma once

#include <memory>

#include "wavefunction.h"

/// @brief The SimpleGaussian class is a gaussian product trial wave function.
/// @details The trial wave function is a product of gaussian functions
/// for each particle, and is independent of the other particles.
class SimpleGaussian : public WaveFunction
{
public:
    /// @brief Constructor for the SimpleGaussian class.
    /// @param alpha Variational parameter present in the exponent.
    SimpleGaussian(double alpha);
    /// @brief Evaluate the trial wave function.
    /// @param particles Vector of particles.
    /// @return The value of the trial wave function.
    double evaluate(std::vector<std::unique_ptr<class Particle>> &particles);
    /// @brief Compute the double derivative of the trial wave function over trial wave function.
    /// @param particles Vector of particles.
    /// @return The local value of Laplasian.
    double computeLocalLaplasian(std::vector<std::unique_ptr<class Particle>> &particles);
    /// @brief Efficiently compute ratio of evaluation of trial wave function.
    /// @param particles_numerator Vector of particles in the numerator.
    /// @param particles_denominator Vector of particles in the denominator.
    double evaluateRatio(std::vector<std::unique_ptr<class Particle>> &particles_numerator,
                         std::vector<std::unique_ptr<class Particle>> &particles_denominator);
};
