#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include <cmath>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega, bool mode) :
        Hamiltonian(system) {
    assert(omega > 0);
	m_mode = mode;
    m_omega  = omega;
}

double HarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {
    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     */

	double r2 = 0;
	for( int dim=0; dim < m_system->getNumberOfDimensions(); dim++ )
	{
		r2 = std::pow(particles[0]->getPosition()[dim], 2);
	}

	double potentialEnergy	= 0.5*m_omega*r2;
	double kineticEnergy = 0;

	//Mulig opptimalisering: computeDoubleDerivative regner ut r2, men vi har allerede regnet 
	//den ut her
	if( m_mode ){
		kineticEnergy	=
			-0.5*m_system->getWaveFunction()->computeDoubleDerivative(particles);
	} else
	{
		kineticEnergy	= numeric();// -.5*numeric();
	}

    return kineticEnergy + potentialEnergy;
}

