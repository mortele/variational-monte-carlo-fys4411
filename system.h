#pragma once

class System {
public:
    System();
    bool metropolisStep             ();
    void runMetropolisSteps         (int numberOfMetropolisSteps);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setStepLength              (double stepLength);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setInitialState            (class InitialState* initialState);
    class WaveFunction* getWaveFunction() { return this->waveFunction; }
    class Hamiltonian*  getHamiltonian()  { return this->hamiltonian; }
    class Particle*     getParticles()    { return this->particles; }
    class Sampler*      getSampler()      { return this->sampler; }
    int getNumberOfParticles()            { return this->numberOfParticles; }
    int getNumberOfDimensions()           { return this->numberOfDimensions; }
    int getNumberOfMetropolisSteps()      { return this->numberOfMetropolisSteps; }
    double getEquilibrationFraction()     { return this->equilibrationFraction; }

private:
    int                     numberOfParticles;
    int                     numberOfDimensions;
    int                     numberOfMetropolisSteps;
    double                  equilibrationFraction = 0.0;
    double                  stepLength = 0.1;
    class WaveFunction*     waveFunction;
    class Hamiltonian*      hamiltonian;
    class Particle*         particles;
    class InitialState*     initialState;
    class Sampler*          sampler;
};

