#include <memory>
#include "wavefunction.h"


WaveFunction::WaveFunction(std::shared_ptr<System> system) {
    m_system = system;
}
