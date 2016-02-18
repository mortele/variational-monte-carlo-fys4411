#include <iostream>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/gaussian4.h"
#include "WaveFunctions/multiparticleho.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
#include "examples.h"

using namespace std;

int main() {
    //System* system = NonInteractingHO(3,20); // (dimensions, particles).
    System* system = InteractingHO(3,20); // (dimensions, particles).
    return 0;
}

