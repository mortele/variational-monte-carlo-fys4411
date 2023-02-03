#pragma once

#include <vector>
#include <memory>

class MonteCarlo {
public:
    MonteCarlo(std::unique_ptr<class Random> rng);
    virtual ~MonteCarlo() = default;

    virtual bool step(
            double stepLength,
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles) = 0;

protected: // originally private?
    std::unique_ptr<class Random> m_rng;
};
