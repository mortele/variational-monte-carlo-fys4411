#include <memory>
#include "hamiltonian.h"
#include "../system.h"

Hamiltonian::Hamiltonian(std::shared_ptr<class System> system) {
    m_system = system;
}

