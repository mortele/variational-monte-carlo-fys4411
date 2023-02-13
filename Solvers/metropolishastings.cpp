#include <memory>
#include <vector>

#include "metropolishastings.h"
#include "WaveFunctions/wavefunction.h"
#include "particle.h"
#include "Math/random.h"

#include <iostream>

MetropolisHastings::MetropolisHastings(std::unique_ptr<class Random> rng)
    : MonteCarlo(std::move(rng))
{
}


bool MetropolisHastings::step(
        double stepLength,
        class WaveFunction& waveFunction,
        std::vector<std::unique_ptr<class Particle>>& particles)
{
    /*
    Here the step related to Metroplis-Hasting algo should be performed, using importance sampling.
    */
   
    static std::vector<double> f(3); // Is static benifical here?

    waveFunction.quantumForce(*particles.at(0), f);

    return true;
}
