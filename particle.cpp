#include "particle.h"
#include <iostream>
#include <cassert>

Particle::Particle() {
}

void Particle::setPosition(const std::vector<double> &position) {
    m_position = position;
}

void Particle::adjustPosition(double change, int dimension) {
    m_position.at(dimension) += change;
}

void Particle::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}
