#pragma once
#include <vector>
#include <chrono>
#include "Math/random.h"

class System {
public:
    System();
    System(int seed);
    bool metropolisStep                     (int particle);
    void runMetropolisSteps                 (int numberOfMetropolisSteps, bool saveData, std::string output);
    bool importanceSamplingStep             (int particle);
    void runImportanceSamplingSteps         (int numberOfMetropolisSteps, bool saveData, std::string output);
    std::vector<double> quantumForce        (int particle);
    void setNumberOfParticles               (int numberOfParticles);
    void setNumberOfDimensions              (int numberOfDimensions);
    void setStepLength                      (double stepLength);
    void setEquilibrationFraction           (double equilibrationFraction);
    void setHamiltonian                     (class Hamiltonian* hamiltonian);
    void setWaveFunction                    (class WaveFunction* waveFunction);
    void setInitialState                    (class InitialState* initialState);
    double GreensFunctionRatio              (std::vector<double> y, std::vector<double> x, 
                                             double dt, std::vector<double> qForceOld, 
                                             std::vector<double> qForceNew);
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    std::vector<class Particle*>    getParticles()      { return m_particles; }
    class Random*                   getRandomEngine()   { return m_random; }
    std::vector<double>             getSdRes()          { return m_sdRes; }
    int getNumberOfParticles()          { return m_numberOfParticles; }
    int getNumberOfDimensions()         { return m_numberOfDimensions; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    double getEquilibrationFraction()   { return m_equilibrationFraction; }
    // std::chrono::time_point<std::chrono::system_clock> getTimeStart()   { return m_time_start;}
    double getElapsedTime()                      { return m_elapsed_time;}

private:
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength = 0.1;
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();
    std::vector<double>             m_sdRes = std::vector<double> {};
    class Random*                   m_random = nullptr;
    std::chrono::time_point<std::chrono::system_clock> m_time_start;
    double m_elapsed_time;
};

