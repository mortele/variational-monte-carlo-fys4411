#include <memory>
#include <iostream>
#include <cassert>

#include "../include/initialstate.h"
#include "../include/particle.h"
#include "../include/random.h"

std::vector<std::unique_ptr<Particle>> setupRandomUniformInitialState(
    double stepLength,
    size_t numberOfDimensions,
    size_t numberOfParticles,
    Random &rng)
{
    assert(numberOfDimensions > 0 && numberOfParticles > 0);

    auto particles = std::vector<std::unique_ptr<Particle>>();

    for (size_t i = 0; i < numberOfParticles; i++)
    {
        std::vector<double> position = std::vector<double>();
        for (size_t j = 0; j < numberOfDimensions; j++)
        {
            //uniformly distributed random number between -stepLength/2 and stepLength/2
            double pos = -stepLength / 2 + rng.nextDouble() * stepLength;
            position.push_back(pos);
        }
        particles.push_back(std::make_unique<Particle>(position));
    }

    return particles;
}
