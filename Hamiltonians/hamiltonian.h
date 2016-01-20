#pragma once

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    void setWaveFunction(class WaveFunction* waveFunction);
    virtual double computeLocalEnergy(class Particle* particles) = 0;

protected:
    class System* m_system = nullptr;

private:
    class WaveFunction* m_waveFunction = nullptr;
};

