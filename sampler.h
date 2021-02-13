#include <string>
#pragma once

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void getOutput();
    void printOutputToTerminal();
    void printOutputToFile();
    void computeAverages();
    double getEnergy()          { return m_energy; }

private:
    int                 m_numberOfMetropolisSteps = 0;
    int                 m_numberOfSampledMetropolisSteps = 0;
    int                 m_stepNumber = 0;
    double              m_energy = 0;
    double              m_cumulativeEnergy = 0;
    std::string         m_output = "";
    std::vector<double> m_energies = std::vector<double>();
    std::vector<double> m_positions = std::vector<double>();
    class System* m_system = nullptr;
};
