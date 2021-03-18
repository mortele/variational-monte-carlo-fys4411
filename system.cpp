#include "system.h"
#include <cassert>
#include <cmath>
#include <chrono>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
using namespace std::chrono; 


System::System() {
    m_random = new Random();
	m_time_start = high_resolution_clock::now();
}

System::System(int seed) {
    m_random = new Random(seed);
	m_time_start = high_resolution_clock::now();
}

std::vector<double> System::quantumForce(int particle) {
	// double waveFun = m_waveFunction->evaluate(m_particles, particle);
	std::vector<double> qForce(m_numberOfParticles);
	std::vector<double> pos = m_particles[particle]->getPosition();
	for (int dim = 0; dim < m_numberOfParticles; dim++)
	{
		qForce[dim] = - 4 * m_waveFunction->getParameters()[0] * pos[dim];
	}
	return qForce;
}

double System::GreensFunctionRatio(std::vector<double> posNew, std::vector<double> posOld,
												double dt, std::vector<double> qForceOld, 
												std::vector<double> qForceNew) {
	/* Function that calculates G(x, y, dt)/G(y, x, dt), where G is
	 * the Green's function, x is the old postion and y is the new position.
	 */
	std::vector<double> gFuncOld(m_numberOfParticles);
	std::vector<double> gFuncNew(m_numberOfParticles);
	double gFunc;
	double gOld = 0;
	double gNew = 0;

	for (int dim = 0; dim < m_numberOfParticles; dim++)
	{
		gFuncOld[dim] = posNew[dim] - posOld[dim] - .5*dt*qForceOld[dim];
		gOld -= gFuncOld[dim] * gFuncOld[dim];
		// gFuncOld[dim] = -gFuncOld[dim]*gFuncOld[dim] / (2*dt);

		gFuncNew[dim] = posOld[dim] - posNew[dim] - .5*dt*qForceNew[dim];
		gNew -= gFuncNew[dim] * gFuncNew[dim];
		// gFuncNew[dim] = - gFuncNew[dim]*gFuncNew[dim] / (2*dt);
		
		// gFunc[dim] = exp(gFuncOld[dim] - gFuncNew[dim]);
	}
	gFunc = exp( (gNew - gOld) * 0.5 / dt );
	return gFunc;
}


bool System::metropolisStep(int particle) {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
	 * double step = m_stepLength*(Random::nextDouble() - .5);
     */

	std::vector<double> step(m_numberOfDimensions);
	double wfold = m_waveFunction->evaluate(m_particles);

	for( int dim = 0; dim < m_numberOfDimensions; dim++ )
	{
		step[dim] = m_stepLength*(m_random->nextDouble() - .5);
		m_particles[particle]->adjustPosition(step[dim], dim);
	}
	double wfnew = m_waveFunction->evaluate(m_particles);

	double ratio = wfnew*wfnew/(wfold*wfold);
	if( m_random->nextDouble() <= ratio )//std::exp(2*(wfnew - wfold)) )
	{
		return true;
	}
	else
	{
		for( int dim = 0; dim < m_numberOfDimensions; dim++ )
		{
			m_particles[particle]->adjustPosition(-step[dim], dim);
		}
		return false;
	}
} //!metropolisStep 

void System::runMetropolisSteps(int numberOfMetropolisSteps, bool saveData) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    for (int i=0; i < numberOfMetropolisSteps; i++) {
		for( int particle=0; particle < m_numberOfParticles; particle++ )
		{
			bool acceptedStep = metropolisStep(particle);

			/* Here you should sample the energy (and maybe other things using
			* the m_sampler instance of the Sampler class. Make sure, though,
			* to only begin sampling after you have let the system equilibrate
			* for a while. You may handle this using the fraction of steps which
			* are equilibration steps; m_equilibrationFraction.
			*/
			if( i > numberOfMetropolisSteps*m_equilibrationFraction )
			{
				m_sampler->sample(acceptedStep);
			}
		}
    }
    m_sampler->computeAverages();
	m_sdRes = {m_sampler->getEnergy(), m_sampler->getDeltaPsi(), m_sampler->getDerivativePsiE()};

	//gets end time
	time_point<system_clock> time_end = high_resolution_clock::now();
	m_elapsed_time = duration_cast<nanoseconds> (time_end - m_time_start).count() / 1e9;
	
	if (saveData)
	{
		m_sampler->getOutput();
		m_sampler->printOutputToTerminal();
		m_sampler->printOutputToFile();
	}
}

bool System::importanceSamplingStep(int particle) {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
	 * double step = m_stepLength*(Random::nextDouble() - .5);
     */

	std::vector<double> step(m_numberOfDimensions);
	double wfold = m_waveFunction->evaluate(m_particles, particle);

	std::vector<double> qForceOld = quantumForce(particle);
	std::vector<double> posOld = m_particles[particle]->getPosition();
	double dt = m_waveFunction->getParameters()[1];
	double xi = m_random->nextGaussian(0, 1);
	for( int dim = 0; dim < m_numberOfDimensions; dim++ )
	{
		step[dim] = .5*qForceOld[dim]*dt + xi*sqrt(dt);
		m_particles[particle]->adjustPosition(step[dim], dim);
	}
	double wfnew = m_waveFunction->evaluate(m_particles, particle);
	std::vector<double> qForceNew = quantumForce(particle);

	double ratio = wfnew*wfnew/(wfold*wfold);
	ratio = ratio*GreensFunctionRatio(m_particles[particle]->getPosition(), posOld, dt, qForceOld, qForceNew);
	if( m_random->nextDouble() <= ratio )//std::exp(2*(wfnew - wfold)) )
	{
		return true;
	}
	else
	{
		for( int dim = 0; dim < m_numberOfDimensions; dim++ )
		{
			m_particles[particle]->adjustPosition(-step[dim], dim);
		}
		return false;
	}
} //!metropolisStep 

void System::runImportanceSamplingSteps(int numberOfMetropolisSteps, bool saveData) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    for (int i=0; i < numberOfMetropolisSteps; i++) {
		for( int particle=0; particle < m_numberOfParticles; particle++ )
		{
			bool acceptedStep = importanceSamplingStep(particle);

			/* Here you should sample the energy (and maybe other things using
			* the m_sampler instance of the Sampler class. Make sure, though,
			* to only begin sampling after you have let the system equilibrate
			* for a while. You may handle this using the fraction of steps which
			* are equilibration steps; m_equilibrationFraction.
			*/
			if( i > numberOfMetropolisSteps*m_equilibrationFraction )
			{
				m_sampler->sample(acceptedStep);
			}
		}
    }
    m_sampler->computeAverages();
	m_sdRes = {m_sampler->getEnergy(), m_sampler->getDeltaPsi(), m_sampler->getDerivativePsiE()};

	//gets end time
	time_point<system_clock> time_end = high_resolution_clock::now();
	m_elapsed_time = duration_cast<nanoseconds> (time_end - m_time_start).count() / 1e9;

    if (saveData)
	{
		m_sampler->getOutput();
		m_sampler->printOutputToTerminal();
		m_sampler->printOutputToFile();
	}
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}


