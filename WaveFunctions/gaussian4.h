#pragma once
#include "wavefunction.h"


class Gaussian4 : public WaveFunction {
public:
    Gaussian4(class System* system, double alpha);
    double evaluate(Particle *particles);
    double computeKineticEnergy(Particle *particles);
};
