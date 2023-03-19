#pragma once

#include <memory>
#include <vector>

class System
{
public:
    System(std::unique_ptr<class Hamiltonian> hamiltonian,
           std::unique_ptr<class WaveFunction> waveFunction,
           std::unique_ptr<class MonteCarlo> solver,
           std::vector<std::unique_ptr<class Particle>> particles);

    unsigned int runEquilibrationSteps(double stepLength,
                                       unsigned int numberOfEquilibrationSteps);

    std::unique_ptr<class Sampler> runMetropolisSteps(
        double stepLength, unsigned int numberOfMetropolisSteps);

    std::unique_ptr<class Sampler> optimizeMetropolis(
        System &system,
        double stepLength,
        unsigned int numberOfMetropolisSteps,
        unsigned int numberOfEquilibrationSteps,
        double epsilon,
        double learningRate);

    double computeLocalEnergy();
    const std::vector<double> &getWaveFunctionParameters();

    void setWaveFunction(std::unique_ptr<class WaveFunction> waveFunction);
    void setSolver(std::unique_ptr<class MonteCarlo> solver);
    double computeParamDerivative(int paramIndex);

private:
    unsigned int m_numberOfParticles = 0;
    unsigned int m_numberOfDimensions = 0;

    std::unique_ptr<class Hamiltonian> m_hamiltonian;
    std::unique_ptr<class WaveFunction> m_waveFunction;
    std::unique_ptr<class MonteCarlo> m_solver;
    std::vector<std::unique_ptr<class Particle>> m_particles;
};
