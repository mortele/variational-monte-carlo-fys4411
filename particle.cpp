#include "particle.h"

Particle::Particle() {
}

void Particle::setPosition(double* position) {
    this->position = position;
}

void Particle::adjustPosition(double change, int dimension) {
    this->position[dimension] += change;
}

void Particle::setNumberOfDimensions(int numberOfDimensions) {
    this->numberOfDimensions = numberOfDimensions;
}
