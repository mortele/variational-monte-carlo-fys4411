#pragma once


class WaveFunction {
public:
    WaveFunction(class System* system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    double* getParameters()         { return m_parameters; }
    virtual double evaluate(class Particle* particles) = 0;
    virtual double computeDoubleDerivative(class Particle* particles) = 0;

protected:
    int     m_numberOfParameters = 0;
    double* m_parameters = nullptr;
    class System* m_system = nullptr;
};

