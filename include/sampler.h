#pragma once
#include <memory>

class Sampler
{
public:
    Sampler(
        size_t numberOfParticles,
        size_t numberOfDimensions,
        double stepLength,
        size_t numberOfMetropolisSteps);

    void sample(bool acceptedStep, class System *system);
    void printOutputToTerminal(class System &system);
    void computeAverages();
    double getEnergy() { return m_energy; }

private:
    size_t m_stepNumber = 0;
    size_t m_numberOfMetropolisSteps = 0;
    size_t m_numberOfParticles = 0;
    size_t m_numberOfDimensions = 0;
    size_t m_numberOfAcceptedSteps = 0;
    double m_energy = 0;
    double m_cumulativeEnergy = 0;
    double m_stepLength = 0;
};
