#pragma once

class InitialState {
public:
    InitialState(class System* system);
    virtual void setupInitialState() = 0;
    class Particle* getParticles() { return m_particles; }

protected:
    class System* m_system;
    class Particle* m_particles;
    int m_numberOfDimensions;
    int m_numberOfParticles;
};

