#pragma once
#include <memory>
#include "initialstate.h"
#include "../system.h"

class RandomUniform : public InitialState {
public:
    RandomUniform(std::shared_ptr<class System> system, unsigned int numberOfDimensions, unsigned int numberOfParticles);
    void setupInitialState();
};

