#include "system.h"

#include <cassert>
#include <iostream>
#include <memory>

#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Solvers/montecarlo.h"
#include "WaveFunctions/wavefunction.h"
#include "particle.h"
#include "sampler.h"

System::System(std::unique_ptr<class Hamiltonian> hamiltonian,
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
    double stepLength, unsigned int numberOfEquilibrationSteps)
{
  unsigned int acceptedSteps = 0;

  // std::cout << m_hamiltonian->computeLocalEnergy(*m_waveFunction,
  // m_particles) << "\n";

  for (unsigned int i = 0; i < numberOfEquilibrationSteps; i++)
  {
    acceptedSteps += m_solver->step(stepLength, *m_waveFunction, m_particles);
  }

  return acceptedSteps;
}

std::unique_ptr<class Sampler> System::runMetropolisSteps(
    double stepLength, unsigned int numberOfMetropolisSteps)
{
  auto sampler =
      std::make_unique<Sampler>(m_numberOfParticles, m_numberOfDimensions,
                                stepLength, numberOfMetropolisSteps);

  for (unsigned int i = 0; i < numberOfMetropolisSteps; i++)
  {
    /* Call solver method to do a single Monte-Carlo step.
     */
    bool acceptedStep =
        m_solver->step(stepLength, *m_waveFunction, m_particles);

    /* Here you should sample the energy (and maybe other things) using the
     * sampler instance of the Sampler class.
     */
    // compute local energy
    sampler->sample(acceptedStep, this);
  }
  sampler->computeAverages();

  return sampler;
}

std::unique_ptr<class Sampler> System::optimizeMetropolis(
    System &system, double stepLength, unsigned int numberOfMetropolisSteps,
    double epsilon, double learningRate)
{
  auto sampler =
      std::make_unique<Sampler>(m_numberOfParticles, m_numberOfDimensions,
                                stepLength, numberOfMetropolisSteps);

  double gradient = 1;
  while (abs(gradient) > epsilon)
  {
    gradient = 0;
    /*Notice that the positions are reset to what they were at the initial state
    but the parameters of the wave function should be what they were at the END
    of last epoch*/

    // reset position and quantum force
    // (i think the quantum force is reset automatically to 0 and the beggining
    // of the metroplis step)
    for (unsigned int i = 0; i < m_numberOfParticles; i++)
    {
      m_particles[i]->resetPosition();
      // m_particles[i]->resetQuantumForce();
    }

    // (re)set the sampler cumulative values by calling the constructor
    sampler =
        std::make_unique<Sampler>(m_numberOfParticles, m_numberOfDimensions,
                                  stepLength, numberOfMetropolisSteps);

    int n_params = m_waveFunction->getNumberOfParameters();

    // call run metropolis steps
    sampler = runMetropolisSteps(stepLength, numberOfMetropolisSteps);

    std::vector<double> m_energyDerivative = sampler->getEnergyDerivative();

    // update parameters
    std::vector<double> parameters = getWaveFunctionParameters();
    for (int i = 0; i < n_params; i++)
    {
      parameters[i] -= learningRate * m_energyDerivative[i]; // gradient descent but this gradient is wrong
      std::cout << "parameters post update: " << parameters[0] << "\n";
      std::cout << "m_energyDerivative: " << m_energyDerivative[i] << "\n";
      gradient += m_energyDerivative[i];
    }
    // set new wave function parameters
    m_waveFunction->setParameters(parameters);
  }
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

double System::computeParamDerivative(int paramIndex)
{
  // Helper function
  return m_waveFunction->computeParamDerivative(m_particles, paramIndex);
}