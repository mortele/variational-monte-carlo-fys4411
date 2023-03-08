#include "particle.h"
#include <cassert>

Particle::Particle(const std::vector<double> &position)
{
    m_numberOfDimensions = position.size();
    m_position = position;
    m_initialPosition = position; // Save initial position. Notice that this is a copy, not a reference.
}

void Particle::adjustPosition(double change, unsigned int dimension)
{
    m_position.at(dimension) += change;
}
void Particle::setPosition(double new_position, unsigned int dimension)
{
    m_position.at(dimension) = new_position;
}
void Particle::resetPosition()
{
    m_position = m_initialPosition;
}