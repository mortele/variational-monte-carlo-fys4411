#pragma once
#include <vector>

class Particle
{
public:
    /// @brief Constructor for the Particle class.
    /// @param position Vector of positions of the particles.
    Particle(const std::vector<double> &position);
    /// @brief Adjust the position of the particle.
    /// @param change The change in position.
    /// @param dimension The dimension to change.
    void adjustPosition(double change, unsigned int dimension);
    /// @brief Get the position of the particle.
    /// @return The position of the particle.
    std::vector<double> getPosition() { return m_position; }
    /// @brief Get the number of dimensions.
    /// @return Number of dimensions 
    unsigned int getNumberOfDimensions() { return m_numberOfDimensions; }

private:
    unsigned int m_numberOfDimensions = 0;
    std::vector<double> m_position = std::vector<double>();
};
