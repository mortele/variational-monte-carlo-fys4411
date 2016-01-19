#pragma once

class InitialState {
public:
    InitialState(class System* system);
    virtual void setupInitialState() = 0;
    class Particle* getParticles() { return this->particles; }

protected:
    class System* system;
    class Particle* particles;
    int numberOfDimensions;
    int numberOfParticles;
};

