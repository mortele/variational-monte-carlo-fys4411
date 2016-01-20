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
    bool    m_exactGroundStateEnergyKnown = 0;
    double  m_exactEnergy = 0;

    double computeKineticEnergy(class Particle* particles);

private:
    class WaveFunction* m_waveFunction = nullptr;
};

