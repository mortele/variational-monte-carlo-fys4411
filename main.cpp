#include <iostream>
#include <fstream>
#include <chrono>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"

#if defined(_WIN32)
    #include <windows.h>
#elif defined(__linux__)
    #include <filesystem>
#endif
// #include <windows.h>
// #include <filesystem>

using namespace std;
using namespace std::chrono;


int main() {
    // Seed for the random number generator
    int seed = 2020;

    int numberOfDimensions[]    = {1, 2};
    int numberOfParticles[]     = {1, 2}; //{1,10,100,500};
    int numberOfSteps           = (int) 1e6;
    double omega                = 1.0;              // Oscillator frequency.
    double alpha[]              = {.4};// Variational parameter.
	//double alpha[] = {0.3, 0.34, 0.38, 0.42, 0.46, 0.5, 0.54, 0.58, 0.62, 0.66, 0.7};
    double stepLength           = 2;              // Metropolis step length.
    double equilibration        = 0.1;              // Amount of the total steps used
    // for equilibration.

    //creares a folder for the results
    #if defined(_WIN32)
        // do some cool Windows stuff
        // #include <windows.h>
        cout << "Windows detected\n";
        string res_folder = ".\\results";
        CreateDirectory(res_folder.c_str(), NULL);
    #elif defined(__linux__)
        // do some cool Unix stuff
        cout << "LINUX detected\n";
        // #include <filesystem>
        namespace fs = std::filesystem;
        if (!fs::is_directory("results") || !fs::exists("results")) { // Check if res folder exists
            fs::create_directory("results"); // create res folder
        }
    #else
        // #error Unsupported operating system
        cout << "No supported os detected. To get results saved in a file create a folder named 'results'.\n";
    #endif 

    // clears output file
    ofstream outfile;
    outfile.open ("results/results.csv", ios::out | ios::trunc);
    outfile << 
    "nParticles;nDimensions;nMetropolisSteps;EquilibrationFraction;acceptedSteps;foundEnergy;elapsedTime;nParameters;Parameters(undefinedNumber)\n";
    outfile.close();
    outfile.open ("results/energies.csv", ios::out | ios::trunc);
    outfile.close();

    time_point<system_clock> tot_time_start = high_resolution_clock::now();
    for (int nDim : numberOfDimensions)
    {
        for (int nPar : numberOfParticles)
        {
            for (double nAlpha : alpha)
            {
                System* system = new System(seed);
                system->setHamiltonian              (new HarmonicOscillator(system, omega, false));
                system->setWaveFunction             (new SimpleGaussian(system, nAlpha));
                system->setInitialState             (new RandomUniform(system, nDim, nPar));
                system->setEquilibrationFraction    (equilibration);
                system->setStepLength               (stepLength);
                system->runMetropolisSteps          (numberOfSteps);
            }
        }
    }
    
    time_point<system_clock> tot_time_end = high_resolution_clock::now();
    double tot_elapsed_time = duration_cast<nanoseconds> (tot_time_end - tot_time_start).count() / 1e9;
    cout << endl;
    cout << " Total time elapsed: " << tot_elapsed_time << " seconds\n";
    cout << endl;

    return 0;
}
