#include "hamiltonian.h"
#include "../system.h"

Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

Hamiltonian::numericDifferentiate()
{

	std::vector<class Particle*> particles = m_system->getParticles();

	double wfcur = m_system->getWaveFunction()->evaluate(particles);
	
	for( int particle=0; 
			particle < m_system->getNumberOfParticles(); 
			particle++)
	{
		
	}
}
