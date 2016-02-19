#pragma once
#include "Hamiltonians/hamiltonian.h"

class HeliumAtom : public Hamiltonian {
public:
    HeliumAtom(class System* system);
    HeliumAtom(class System* system, bool interaction);
    double computeLocalEnergy(Particle *particles);

private:
    bool m_interaction = true;
};

