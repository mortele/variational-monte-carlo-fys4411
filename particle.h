#pragma once
#include <vector>

class Particle {
public:
    Particle(const std::vector<double>& position);
    void adjustPosition(double change, unsigned int dimension);
    std::vector<double> &getPosition() { return m_position; }
    unsigned int getNumberOfDimensions() { return m_numberOfDimensions; }

private:
    unsigned int m_numberOfDimensions = 0;
    std::vector<double> m_position = std::vector<double>();
};

