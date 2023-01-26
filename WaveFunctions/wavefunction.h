#pragma once
#include <memory>
#include <vector>


class WaveFunction {
public:
    virtual ~WaveFunction() = default;

    int getNumberOfParameters() { return m_numberOfParameters; }
    const std::vector<double>& getParameters() { return m_parameters; }
    virtual double evaluate(std::vector<std::unique_ptr<class Particle>>& particles) = 0;
    virtual double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles) = 0;

protected:
    int m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
};

