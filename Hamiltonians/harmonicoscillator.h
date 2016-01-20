#pragma once
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega);
    double computeLocalEnergy(std::vector<Particle*> particles);

private:
    double m_omega = 0;
};

