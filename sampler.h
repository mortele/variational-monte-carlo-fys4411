#pragma once

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    double getEnergy()          { return this->energy; }
    double getEnergyVariance()  { return this->energyVariance; }

private:
    int     numberOfMetropolisSteps;
    int     stepNumber;
    double  energy;
    double  energyVariance;
    double  acceptanceRate;
    double  cumulativeEnergy;
    double  cumulativeEnergy2;
    double  cumulativeAcceptanceRate;
    class System* system;
};
