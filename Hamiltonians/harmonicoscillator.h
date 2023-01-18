#pragma once
#include <memory>
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(std::shared_ptr<class System> system, double omega);
    double computeLocalEnergy(std::vector<std::unique_ptr<Particle>> particles);

private:
    double m_omega;
};

