#include "hamiltonian.h"
#include "../system.h"

Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

//double Hamiltonian::numericKinetic()
//{
//
//	std::vector<class Particle*> particles = m_system->getParticles();
//
//	double wfcur = m_system->getWaveFunction()->evaluate(particles);
//	double wfnext = 0; double wfprev = 0;
//	// found in 
//	// ComputationalPhysics2/doc/Programs/LecturePrograms/programs/vmc_atoms.py
//	// lines 91-92:
//	double steplength = 1e-3;
//	double inversesquare = 1./(steplength*steplength);
//
//	double kinetic=0;
//
//	for(	int particle=0; 
//			particle < m_system->getNumberOfParticles(); 
//			particle++)
//	{
//		for(	int dimension=0;
//				dimension < m_system->getNumberOfDimensions();
//				dimension++)
//		{
//			particles[particle]->adjustPosition(steplength, dimension);
//			fwnext = m_system->getWaveFunction()->evaluate(particles);
//			particles[particle]->adjustPosition(-2*steplength, dimension);
//
//			wfprev = m_system->getWaveFunction()->evaluate(particles);
//			particles[particle]->adjustPosition(steplength, dimension);
//
//			kinetic -= fwnext + wfprev - 2*wfcur;
//		}//!dimension
//	}//!particle
//	kinetic = .5*inversesquare*kinetic/wfcur;
//	return kinetic;
//}
