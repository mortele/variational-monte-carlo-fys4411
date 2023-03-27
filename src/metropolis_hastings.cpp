#include <memory>
#include <vector>

#include "../include/metropolis_hastings.h"
#include "../include/wavefunction.h"
#include "../include/particle.h"
#include "../include/random.h"

MetropolisHastings::MetropolisHastings(std::unique_ptr<class Random> rng)
    : MonteCarlo(std::move(rng))
{
}

bool MetropolisHastings::step(
    double stepLength,
    class WaveFunction &waveFunction,
    std::vector<std::unique_ptr<class Particle>> &particles)
{
	double D = 0.5; //Parameter in the Fokker-Plank simulation
	
	double rootOfStepLength = sqrt(stepLength); 
	
    // Choose a random particle
    auto particle_index = m_rng->nextInt(particles.size() - 1);
    auto particle = *particles[particle_index];

	//Calculate quantum force for particle before change
	auto quantumForce = waveFunction.computeQuantumForce(particles, particle_index);

    // Save old particle
    auto new_particle = particle;

    // Propose a new position, using Langevin equation (think of it as motion due to force plus some random noise).
    for (size_t pos_index = 0; pos_index < particle.getNumberOfDimensions(); ++pos_index)
        new_particle.adjustPosition(0.5*quantumForce[pos_index]*stepLength + m_rng->nextGaussian(0,1)*rootOfStepLength, pos_index);

	//Calculate quantum force for particle after change
	auto quantumForceNew = waveFunction.computeQuantumForce(particles, particle_index);

	//Calculate the Greens function
	double green = 0.0;
	for (int j = 0; j < particle.getNumberOfDimensions(); j++){
		green += 0.5*(quantumForce[j] + quantumForceNew[j])*
				(D*stepLength*0.5*(quantumForce[j] - quantumForceNew[j]) - new_particle.getPosition()[j] + particle.getPosition()[j]);
	}
	green = exp(green);
	
    // copy particles array 
    auto new_particles = std::vector<std::unique_ptr<class Particle>>(particles.size());
    for (size_t i = 0; i < particles.size(); ++i)
        new_particles[i] = std::make_unique<Particle>(*particles[i]);

    new_particles[particle_index] = std::make_unique<Particle>(new_particle);

    // Calculate the ratio of the new and old wave functions
    auto ratio = waveFunction.evaluateRatio(new_particles, particles);

    // Accept or reject the move, now using Greens function for a more informed decision, making it now a Metropolis-Hastings algorithm.
    if (m_rng->nextDouble() < ratio * ratio * green)
    {
        // Accept the move
        particles[particle_index] = std::make_unique<Particle>(new_particle);
        return true;
    }

    return false;
}
