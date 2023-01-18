#pragma once
#include <vector>

class Particle {
public:
    Particle(const std::vector<double>& position, unsigned int numberOfDimensions);
    void adjustPosition(double change, unsigned int dimension);
    std::vector<double> getPosition() { return m_position; }

private:
    unsigned int m_numberOfDimensions = 0;
    std::vector<double> m_position = std::vector<double>();
};

