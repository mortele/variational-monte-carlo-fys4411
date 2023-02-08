#include <memory>
#include <iostream>
#include <cassert>

#include "initialstate.h"
#include "../particle.h"
#include "Math/random.h"


std::vector<std::unique_ptr<Particle>> setupRandomUniformInitialState(
            double stepLength,
            unsigned int numberOfDimensions,
            unsigned int numberOfParticles,
            Random& rng
        )
{
    assert(numberOfDimensions > 0 && numberOfParticles > 0);

    auto particles = std::vector<std::unique_ptr<Particle>>();

    for (unsigned int i=0; i < numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();

        for (unsigned int j=0; j < numberOfDimensions; j++) {
            /* This is where you should actually place the particles in
             * some positions, according to some rule. Since this function is
             * called random uniform, they should be placed randomly according
             * to a uniform distribution here.
             *
             * Note: For now, the particles are simply placed in positions
             * according to their index in the particles list (this is
             * NOT a good idea).
             */
            double q = rng.nextDouble()-0.5;
            position.push_back(q);
        }

        particles.push_back(std::make_unique<Particle>(position));
    }

    return particles;
}
