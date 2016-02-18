#pragma once
#include "hamiltonian.h"

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega);
    virtual double computeLocalEnergy(Particle* particles);

protected:
    double m_omega  = 0;
    double m_omega2 = 0;
};

