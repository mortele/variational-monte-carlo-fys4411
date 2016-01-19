#pragma once


class WaveFunction {
public:
    WaveFunction(class System* system);
    int     getNumberOfParameters() { return this->numberOfParameters; }
    double* getParameters()         { return this->parameters; }
    virtual double evaluate(class Particle* particles) = 0;
    virtual double computeDoubleDerivative(class Particle* particles) = 0;

protected:
    int     numberOfParameters;
    double* parameters;
    class System* system;
};

