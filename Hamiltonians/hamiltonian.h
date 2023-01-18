#pragma once
#include <memory>
#include <vector>

class Hamiltonian {
public:
    Hamiltonian(std::shared_ptr<class System> system);
    virtual ~Hamiltonian() = default;

    virtual double computeLocalEnergy(std::vector<std::unique_ptr<class Particle>> particles) = 0;

protected:
    std::shared_ptr<class System> m_system;
};

