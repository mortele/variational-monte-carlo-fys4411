#include "../include/particle.h"
#include <cassert>

Particle::Particle(const std::vector<double> &position)
{
    m_numberOfDimensions = position.size();
    m_position = position;
}

void Particle::adjustPosition(double change, size_t dimension)
{
    m_position.at(dimension) += change;
}
