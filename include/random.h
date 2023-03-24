#pragma once

#include <random>

/* This is a proposal for a convenience random number generator. However, feel
 * free to use the standard library, or any other libabry, to create your own
 * generators. Example usage of this class is:
 *
 *  auto rng = std::make_unique<Random>(2020); // Create a unique rng instance
 *  int foo = rng->nextInt(5, 10); // Draw random uniform integers in [5, 10]
 *  double bar = rng->nextDouble(); // Draw random uniform doubles in [0, 1)
 */

/// @brief Convenience class for random number generation.
class Random
{

private:
    /// @brief The random number generator engine.
    std::mt19937_64 m_engine;

public:
    /// @brief Constructor for random number generator.
    Random()
    {
        std::random_device rd;
        m_engine = std::mt19937_64(rd());
    }

    /// @brief Constructor for random number generator.
    Random(int seed)
    {
        m_engine = std::mt19937_64(seed);
    }

    /// @brief Produces uniformly distributed random integers in the closed interval [lowerLimit, upperLimit].
    /// @param lowerLimit lower limit of the interval
    /// @param upperLimit upper limit of the interval
    /// @return random integer in the interval [lowerLimit, upperLimit]
    int nextInt(int lowerLimit, int upperLimit)
    {
        std::uniform_int_distribution<int> dist(lowerLimit, upperLimit);
        return dist(m_engine);
    }

    /// @brief Produces uniformly distributed random integers in the closed interval [0, upperLimit].
    /// @param upperLimit upper limit of the interval
    /// @return random integer in the interval [0, upperLimit]
    int nextInt(int upperLimit)
    {
        std::uniform_int_distribution<int> dist(0, upperLimit);
        return dist(m_engine);
    }

    /// @brief Produces uniformly distributed random floating-point values in the half-open interval [0, 1).
    /// @return random double in the interval [0, 1)
    double nextDouble()
    {
        std::uniform_real_distribution<double> dist(0, 1);
        return dist(m_engine);
    }

    /// @brief Produces gaussian distributed random floating-point values with mean ``mean`` and standard deviation ``standardDeviation``.
    /// @param mean mean of the distribution
    /// @param standardDeviation standard deviation of the distribution
    /// @return random double with mean ``mean`` and standard deviation ``standardDeviation``
    double nextGaussian(
        double mean,
        double standardDeviation)
    {
        std::normal_distribution<double> dist(mean, standardDeviation);
        return dist(m_engine);
    }
};
