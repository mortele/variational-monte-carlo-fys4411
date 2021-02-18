#pragma once

#include <random>

/* This is a proposal for a convenience random number generator. However, feel
 * free to use the standard library, or any other libabry, to create your own
 * generators. Example usage of this class is:
 *
 *  class Random* rng = new Random(2020); // Create rng instance
 *  int foo = rng->nextInt(5, 10); // Draw random uniform integers in [5, 10]
 *  double bar = rng->nextDouble(); // Draw random uniform doubles in [0, 1)
 */


class Random {

private:
    std::mt19937_64 m_engine;

public:
    Random() {
        std::random_device rd;
        m_engine = std::mt19937_64(rd());
    }

    Random(int seed) {
        m_engine = std::mt19937_64(seed);
    }

    int nextInt(const int &lowerLimit, const int &upperLimit) {
        // Produces uniformly distributed random integers in the closed interval
        // [lowerLimit, upperLimit].

        std::uniform_int_distribution<int> dist(lowerLimit, upperLimit);
        return dist(m_engine);
    }

    int nextInt(const int &upperLimit) {
        // Produces uniformly distributed random integers in the closed interval
        // [0, upperLimit].

        std::uniform_int_distribution<int> dist(0, upperLimit);
        return dist(m_engine);
    }

    double nextDouble() {
        // Produces uniformly distributed random floating-point values in the
        // half-open interval [0, 1).

        std::uniform_real_distribution<double> dist(0, 1);
        return dist(m_engine);
    }
    double nextGaussian(
        const double &mean,
        const double &standardDeviation
    ) {
        // Produces normal distributed random floating-point values with mean
        // ``mean`` and standard deviation ``standardDeviation``.

        std::normal_distribution<double> dist(mean, standardDeviation);
        return dist(m_engine);
    }
};

