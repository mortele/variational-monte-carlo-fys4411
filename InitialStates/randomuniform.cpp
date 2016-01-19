#include "randomuniform.h"
#include "Math/random.h"
#include "../particle.h"
#include "../system.h"


RandomUniform::RandomUniform(System*    system,
                             int        numberOfDimensions,
                             int        numberOfParticles)  :
        InitialState(system) {
    this->numberOfDimensions = numberOfDimensions;
    this->numberOfParticles  = numberOfParticles;
    this->system->setNumberOfDimensions(numberOfDimensions);
    this->system->setNumberOfParticles(numberOfParticles);
    this->setupInitialState();
}

void RandomUniform::setupInitialState() {
    this->particles = new Particle[this->numberOfParticles];

    for (int i=0; i < this->numberOfParticles; i++) {
        double* position = new double[this->numberOfDimensions];
        for (int j=0; j < this->numberOfDimensions; j++) {
            position[j] = Random::nextDouble()*2-1;
        }
        this->particles[i].setNumberOfDimensions(this->numberOfDimensions);
        this->particles[i].setPosition(position);
    }
}
