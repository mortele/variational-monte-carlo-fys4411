#pragma once
#include <memory>
#include <vector>

#include "particle.h"

/// @brief The WaveFunction class is the base class for all trial wave functions.
/// @details It contains the number of parameters and the parameters themselves.
class WaveFunction
{
public:
    /// @brief Destructor for the WaveFunction class.
    virtual ~WaveFunction() = default;

    /// @brief Get the number of parameters.
    /// @return The number of parameters.
    size_t getNumberOfParameters() const { return m_numberOfParameters; }
    /// @brief Get the const reference to parameters.
    /// @return The parameter list.
    const std::vector<double> &getParameters() const { return m_parameters; }
    /// @brief Evaluate the trial wave function.
    /// @param particles Vector of particles.
    /// @return The value of the trial wave function.
    virtual double evaluate(std::vector<std::unique_ptr<class Particle>> &particles) = 0;
    /// @brief Compute the double derivative of the trial wave function over trial wave function.
    /// @param particles Vector of particles.
    /// @return The local value of Laplasian.
    virtual double computeLocalLaplasian(std::vector<std::unique_ptr<class Particle>> &particles);
    /// @brief Compute ratio of evaluation of trial wave function.
    /// @param particles_numerator Vector of particles in the numerator.
    /// @param particles_denominator Vector of particles in the denominator.
    virtual double evaluateRatio(std::vector<std::unique_ptr<class Particle>> &particles_numerator,
                                 std::vector<std::unique_ptr<class Particle>> &particles_denominator);

    virtual std::vector<double> computeQuantumForce(std::vector<std::unique_ptr<class Particle>> &particles, size_t particle_index);

protected:
    size_t m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
};
