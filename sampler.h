#pragma once
#include <memory>
#include <string>

class Sampler
{
public:
    Sampler(
        unsigned int numberOfParticles,
        unsigned int numberOfDimensions,
        double stepLength,
        unsigned int numberOfMetropolisSteps);

    void sample(bool acceptedStep, class System *system);
    void printOutputToTerminal(class System &system);
    void writeOutToFile(class System &system, std::string filename, double omega, bool analytical);
    void computeAverages();
    double getEnergy()
    {
        return m_energy;
    }

private:
    double m_stepLength = 0;
    unsigned int m_stepNumber = 0;
    unsigned int m_numberOfMetropolisSteps = 0;
    unsigned int m_numberOfParticles = 0;
    unsigned int m_numberOfDimensions = 0;
    unsigned int m_numberOfAcceptedSteps = 0;

    double m_energy = 0;
    double m_energy_variance = 0;
    double m_energy_std = 0;
    double m_acceptRatio = 0;

    double m_cumulativeEnergy = 0;
    double m_cumulativeEnergy2 = 0;
};
