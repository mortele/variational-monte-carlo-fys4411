#pragma once
#include <memory>
#include <vector>
#include "../particle.h"
#include "../system.h"

class InitialState {
public:
    InitialState(std::shared_ptr<class System> system);
    virtual ~InitialState() = default;

    virtual void setupInitialState() = 0;
    std::vector<std::unique_ptr<Particle>> getParticles() { return std::move(m_particles); }

protected:
    std::shared_ptr<class System> m_system = nullptr;
    std::vector<std::unique_ptr<class Particle>> m_particles = std::vector<std::unique_ptr<class Particle>>();
    unsigned int m_numberOfDimensions = 0;
    unsigned int m_numberOfParticles = 0;
};

