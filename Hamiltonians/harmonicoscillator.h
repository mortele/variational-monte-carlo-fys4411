#pragma once
#include "hamiltonian.h"

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega);
    double computeLocalEnergy(Particle* particles);

private:
    double omega;
    double omega2;
};

