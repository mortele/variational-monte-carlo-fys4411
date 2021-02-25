#include "hamiltonian.h"
#include <cmath>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

double Hamiltonian::numeric()
{

	std::vector<class Particle*> particles = m_system->getParticles();

	double wfcur = 0; 
	for (int particle = 0; particle < m_system->getNumberOfParticles(); particle++)
	{
		wfcur += m_system->getWaveFunction()->evaluate(particles, particle);
	}

	double wfnext = 0; double wfprev = 0;
	// found in 
	// ComputationalPhysics2/doc/Programs/LecturePrograms/programs/vmc_atoms.py
	// lines 91-92:
	double h = 1e-3;
	double h2 = 1./(h*h);
	//double h2 = 1e6;

	double deriv=0;

	for(	int particle=0; 
			particle < m_system->getNumberOfParticles(); 
			particle++)
	{
		for(	int dimension=0;
				dimension < m_system->getNumberOfDimensions();
				dimension++)
		{
			particles[particle]->adjustPosition(h, dimension);
			wfnext = m_system->getWaveFunction()->evaluate(particles, particle);
			particles[particle]->adjustPosition(-2*h, dimension);

			wfprev = m_system->getWaveFunction()->evaluate(particles, particle);
			particles[particle]->adjustPosition(h, dimension);

			deriv -= wfnext + wfprev - 2*wfcur;
		}//!dimension
	}//!particle
	return  .5*h2*deriv/wfcur;
} // !numeric
