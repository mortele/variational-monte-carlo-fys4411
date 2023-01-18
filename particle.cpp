#include "particle.h"
#include <cassert>

Particle::Particle(const std::vector<double>& position, unsigned int numberOfDimensions) {
    assert (position.size() == numberOfDimensions);
    m_numberOfDimensions = numberOfDimensions;
    m_position = position;
}

void Particle::adjustPosition(double change, unsigned int dimension) {
    m_position.at(dimension) += change;
}
