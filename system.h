#pragma once
#include <memory>
#include <vector>
#include <Math/random.h>
#include "particle.h"

class System {
public:
    System();
    System(int seed);
    bool metropolisStep();
    void runMetropolisSteps(unsigned int numberOfMetropolisSteps);
    void setNumberOfParticles(unsigned int numberOfParticles);
    void setNumberOfDimensions(unsigned int numberOfDimensions);
    void setStepLength(double stepLength);
    void setEquilibrationFraction(double equilibrationFraction);

    void setHamiltonian(std::unique_ptr<class Hamiltonian> hamiltonian);
    void setWaveFunction(std::unique_ptr<class WaveFunction> waveFunction);
    void setInitialState(std::unique_ptr<class InitialState> initialState);
    // class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    // class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    // class Sampler*                  getSampler()        { return m_sampler; }
    // std::vector<class Particle*>    getParticles()      { return m_particles; }
    // class Random*                   getRandomEngine()   { return m_random; }
    unsigned int getNumberOfParticles() { return m_numberOfParticles; }
    unsigned int getNumberOfDimensions() { return m_numberOfDimensions; }
    unsigned int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    double getEquilibrationFraction() { return m_equilibrationFraction; }

private:
    unsigned int m_numberOfParticles = 0;
    unsigned int m_numberOfDimensions = 0;
    unsigned int m_numberOfMetropolisSteps = 0;
    double m_equilibrationFraction = 0.0;
    double m_stepLength = 0.1;
    std::unique_ptr<class WaveFunction> m_waveFunction;
    std::unique_ptr<class Hamiltonian> m_hamiltonian;
    std::unique_ptr<class InitialState> m_initialState;
    std::unique_ptr<class Sampler> m_sampler;
    std::unique_ptr<Random> m_random;

    std::vector<std::unique_ptr<class Particle>> m_particles = std::vector<std::unique_ptr<Particle>>();
};

