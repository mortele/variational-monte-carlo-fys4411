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
    void writeOutToFile(class System &system, std::string filename, double omega, bool analytical, bool importanceSampling);
    void WriteTimingToFiles(System &system, std::string filename, bool analytical, unsigned int numberOfEquilibrationSteps, double timing);
    void writeGradientSearchToFile(System &system, std::string filename, double alpha_0, int epoch, double alpha, double beta);
    void output(System &system, std::string filename, double omega, bool analytical, bool importanceSampling);
    void computeAverages();
    std::vector<double> getEnergyDerivative();

    double getEnergy()
    {
        return m_energy;
    }

private:
    double m_stepLength = 0;
    unsigned int m_stepNumber = 0;
    unsigned int m_numberOfMetropolisSteps = 0;
    unsigned int m_numberOfEquilibrationSteps = 0;
    unsigned int m_numberOfParticles = 0;
    unsigned int m_numberOfDimensions = 0;
    unsigned int m_numberOfAcceptedSteps = 0;

    double m_energy = 0;
    double m_energy_variance = 0;
    double m_energy_std = 0;
    double m_acceptRatio = 0;

    double m_cumulativeEnergy = 0;
    double m_cumulativeEnergy2 = 0;

    int m_numberOfParams = 1; // this should not be hard coded but we will change it later

    std::vector<double> m_energyDerivative = std::vector<double>(m_numberOfParams, 0);
    std::vector<double> m_cumulativeDerPsiE = std::vector<double>(m_numberOfParams, 0);
    std::vector<double> m_cumulativedeltaPsi = std::vector<double>(m_numberOfParams, 0);
    std::vector<double> m_deltaPsi = std::vector<double>(m_numberOfParams, 0);
    std::vector<double> m_derPsiE = std::vector<double>(m_numberOfParams, 0);
};
