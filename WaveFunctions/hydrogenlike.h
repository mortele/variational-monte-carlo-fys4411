#pragma once
#include "WaveFunctions/wavefunction.h"

class HydrogenLike : public WaveFunction {
public:
    HydrogenLike(class System* system, double alpha, double beta);
    double evaluate(Particle *particles);
    double computeKineticEnergy(Particle *particles);

private:
    int m_numberOfDimensions = 0;
};

