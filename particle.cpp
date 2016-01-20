#include "particle.h"

Particle::Particle() {
}

void Particle::setPosition(double* position) {
    m_position = position;
}

void Particle::adjustPosition(double change, int dimension) {
    m_position[dimension] += change;
}

void Particle::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}
