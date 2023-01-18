#pragma once
#include <memory>

class Sampler {
public:
    Sampler(std::shared_ptr<class System> system);
    void setNumberOfMetropolisSteps(unsigned int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    double getEnergy() { return m_energy; }

private:
    unsigned int m_numberOfMetropolisSteps = 0;
    int m_stepNumber = 0;
    double m_energy = 0;
    double m_cumulativeEnergy = 0;
    std::shared_ptr<class System> m_system;
};
