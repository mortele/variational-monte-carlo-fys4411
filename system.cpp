#include <iostream>
#include <memory>
#include <cassert>

#include "system.h"
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Solvers/montecarlo.h"

#include <iostream>

System::System(
    std::unique_ptr<class Hamiltonian> hamiltonian,
    std::unique_ptr<class WaveFunction> waveFunction,
    std::unique_ptr<class MonteCarlo> solver,
    std::vector<std::unique_ptr<class Particle>> particles)
{
    m_numberOfParticles = particles.size();
    m_numberOfDimensions = particles[0]->getNumberOfDimensions();
    m_hamiltonian = std::move(hamiltonian);
    m_waveFunction = std::move(waveFunction);
    m_solver = std::move(solver);
    m_particles = std::move(particles);
}

unsigned int System::runEquilibrationSteps(
    double stepLength,
    unsigned int numberOfEquilibrationSteps)
{
    unsigned int acceptedSteps = 0;

    // std::cout << m_hamiltonian->computeLocalEnergy(*m_waveFunction, m_particles) << "\n";

    for (unsigned int i = 0; i < numberOfEquilibrationSteps; i++)
    {
        acceptedSteps += m_solver->step(stepLength, *m_waveFunction, m_particles);
    }

    return acceptedSteps;
}

std::unique_ptr<class Sampler> System::runMetropolisSteps(
    double stepLength,
    unsigned int numberOfMetropolisSteps,
    bool gradientDescent,
    double learningRate,
    int epochs)
{
    auto sampler = std::make_unique<Sampler>(
        m_numberOfParticles,
        m_numberOfDimensions,
        stepLength,
        numberOfMetropolisSteps);
    double localEnergy = 0;
    for (int i = 0; i < epochs; i++)
    {
        // Here I think I should reset the position to essentially start over a new metroplis run
        // maybe something like this:
        // for (int i = 0; i < m_numberOfParticles; i++)
        // {
        //     m_particles[i]->resetPosition();
        // }

        for (unsigned int i = 0; i < numberOfMetropolisSteps; i++)
        {
            /* Call solver method to do a single Monte-Carlo step.
             */
            bool acceptedStep = m_solver->step(stepLength, *m_waveFunction, m_particles);

            /* Here you should sample the energy (and maybe other things) using the
             * sampler instance of the Sampler class.
             */
            // compute local energy

            sampler->sample(acceptedStep, this);
        }
        localEnergy = computeLocalEnergy();

        std::cout << "Epoch " << i << "Local energy: " << localEnergy << "\n";

        if (gradientDescent)
        {
            // compute gradient
            std::vector<double> gradient = m_waveFunction->computeDerivative(m_particles);
            // get current parameters
            std::vector<double> parameters = getWaveFunctionParameters();
            // update parameters
            for (int i = 0; i < parameters.size(); i++)
            {
                std::cout << "Parameter " << i << ": " << parameters[i] << "\n";
                std::cout << "Gradient " << i << ": " << gradient[i] << "\n";
                parameters[i] -= learningRate * gradient[i];
            }
            // set new wave function
            m_waveFunction->setParameters(parameters);
        }
    }

    sampler->computeAverages();

    return sampler;
}

double System::computeLocalEnergy()
{
    // Helper function
    return m_hamiltonian->computeLocalEnergy(*m_waveFunction, m_particles);
}

const std::vector<double> &System::getWaveFunctionParameters()
{
    // Helper function
    return m_waveFunction->getParameters();
}

void System::setWaveFunction(std::unique_ptr<class WaveFunction> waveFunction)
{
    m_waveFunction = std::move(waveFunction);
}

void System::setSolver(std::unique_ptr<class MonteCarlo> solver)
{
    m_solver = std::move(solver);
}