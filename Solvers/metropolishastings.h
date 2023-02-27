#pragma once

#include <memory>

#include "montecarlo.h"


class MetropolisHastings : public MonteCarlo {
public:
    MetropolisHastings(std::unique_ptr<class Random> rng, double timeStep, double D);
    bool step(
            double stepLength,
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles);

private:
    double m_timeStep; 
    double m_sqrtTimeStep;
    double m_D;
};