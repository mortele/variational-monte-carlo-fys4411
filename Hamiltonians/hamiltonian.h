#pragma once

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    void setWaveFunction(class WaveFunction* waveFunction);
    virtual double computeLocalEnergy(class Particle* particles) = 0;

protected:
    class System* system;
    double computeKineticEnergy(class Particle* particles);


private:
    class WaveFunction* waveFunction;
};

