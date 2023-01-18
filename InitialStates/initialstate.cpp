#include <memory>
#include "initialstate.h"

InitialState::InitialState(std::shared_ptr<System> system) {
    m_system = system;
}

