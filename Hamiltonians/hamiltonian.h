#pragma once

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    void setWaveFunction(class WaveFunction* waveFunction);
    virtual double computeLocalEnergy(class Particle* particles) = 0;
    double getExactEnergy();
    bool getExactGroundStateEnergyKnown() { return m_exactGroundStateEnergyKnown; }

protected:
    class System* m_system = nullptr;
    class WaveFunction* m_waveFunction = nullptr;
    bool    m_exactGroundStateEnergyKnown = false;
    double  m_exactEnergy = 0;

    virtual double computeKineticEnergy(class Particle* particles);
};

