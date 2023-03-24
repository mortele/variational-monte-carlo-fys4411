#pragma once
#include <vector>

class Particle
{
public:
    Particle(const std::vector<double> &position);
    void adjustPosition(double change, unsigned int dimension);
    void setPosition(double new_position, unsigned int dimension);
    void resetPosition();
    void saveEquilibrationPosition();
    void resetEquilibrationPosition();
    std::vector<double> getPosition() { return m_position; }
    unsigned int getNumberOfDimensions() { return m_numberOfDimensions; }

private:
    unsigned int m_numberOfDimensions = 0;
    std::vector<double> m_position = std::vector<double>();
    std::vector<double> m_initialPosition = std::vector<double>(); // Save initial position to reset in Gradient Descent.
    std::vector<double> m_EquilibrationPosition = std::vector<double>();
};

double particle_r2(Particle &p);
double particle_r2(Particle &p1, Particle &p2);