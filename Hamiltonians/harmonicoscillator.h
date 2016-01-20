#pragma once
#include "hamiltonian.h"

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega);
    double computeLocalEnergy(Particle* particles);

private:
    double m_omega = 0;
};

