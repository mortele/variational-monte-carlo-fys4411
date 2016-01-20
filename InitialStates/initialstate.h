#pragma once

class InitialState {
public:
    InitialState(class System* system);
    virtual void setupInitialState() = 0;
    class Particle* getParticles() { return m_particles; }

protected:
    class System* m_system = nullptr;
    class Particle* m_particles = nullptr;
    int m_numberOfDimensions = 0;
    int m_numberOfParticles = 0;
};

