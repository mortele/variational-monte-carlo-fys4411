#include <memory>
#include <vector>

#include "metropolishastings.h"
#include "WaveFunctions/wavefunction.h"
#include "particle.h"
#include "Math/random.h"

#include <iostream>

MetropolisHastings::MetropolisHastings(std::unique_ptr<class Random> rng, double timeStep, double D)
    : MonteCarlo(std::move(rng))
{
    m_timeStep = timeStep;
    m_sqrtTimeStep = std::sqrt(timeStep);
    m_D = D;
}


bool MetropolisHastings::step(
        double stepLength,
        class WaveFunction& waveFunction,
        std::vector<std::unique_ptr<class Particle>>& particles)
{
    /*
    Here the step related to Metroplis-Hasting algo should be performed, using importance sampling.
    */   
    using namespace std;
    int numberOfParticles = particles.size();
    int numberOfDimensions = particles.at(0)->getNumberOfDimensions();

    std::vector<double> qForceOld(numberOfDimensions);
    std::vector<double> qForceNew(numberOfDimensions);

    int proposed_particle_idx = m_rng->nextInt(0,numberOfParticles-1);
    
    Particle& proposed_particle = *particles.at(proposed_particle_idx);
    Particle old_particle = proposed_particle; 

    double Psi_old = waveFunction.evaluate(particles);
    waveFunction.quantumForce(proposed_particle, qForceOld);

    for(int q = 0; q < numberOfDimensions; q++) {
        proposed_particle.adjustPosition(
            m_D*qForceOld.at(q)*m_timeStep + m_rng->nextGaussian(0.0, 1.0)*m_sqrtTimeStep,
            q
        );
    }

    double Psi_new = waveFunction.evaluate(particles);
    waveFunction.quantumForce(proposed_particle, qForceNew);

    double G_ratio = 0;
    for(int q = 0; q < numberOfDimensions; q++) {
        G_ratio = 0.5*(qForceOld.at(q)+qForceNew.at(q))*(
            old_particle.getPosition().at(q) - proposed_particle.getPosition().at(q)
            + 0.5*m_D*m_timeStep*(qForceOld.at(q) - qForceNew.at(q))
        );
    }
    G_ratio = std::exp(G_ratio);
    
    double w = G_ratio * (Psi_new*Psi_new)/(Psi_old*Psi_old);

    if( w >= m_rng->nextDouble() ) {
        return true;
    }
    else {
        *particles.at(proposed_particle_idx) = old_particle;
        return false;
    }
}