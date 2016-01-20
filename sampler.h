#pragma once

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    double getEnergy()          { return m_energy; }
    double getEnergyVariance()  { return m_energyVariance; }

private:
    int     m_numberOfMetropolisSteps;
    int     m_stepNumber;
    double  m_energy;
    double  m_energyVariance;
    double  m_acceptanceRate;
    double  m_cumulativeEnergy;
    double  m_cumulativeEnergy2;
    double  m_cumulativeAcceptanceRate;
    class System* m_system;
};
