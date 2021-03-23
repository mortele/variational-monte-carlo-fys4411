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
#include <cmath>

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

    int numberOfDimensions[]    = {3};//= {1, 2, 3};
    int numberOfParticles[]     = {1};//= {1, 2, 3}; //{1,10,100,500};
    int numberOfSteps           = (int) std::pow(2, 9);
    double omega                = 1.0;              // Oscillator frequency.
    double alpha[]              = {.4};// Variational parameter.
	//double alpha[] = {0.3, 0.34, 0.38, 0.42, 0.46, 0.5, 0.54, 0.58, 0.62, 0.66, 0.7};
    double stepLength           = 2;              // Metropolis step length.
    double equilibration        = 0.1;              // Amount of the total steps used
    // for equilibration.
    double dt                   = 0.01;

    //creares a folder for the results
    #if defined(_WIN32)
        // do some cool Windows stuff
        cout << "Windows detected\n";
        string res_folder = ".\\results";
        CreateDirectory(res_folder.c_str(), NULL);
    #elif defined(__linux__)
        // do some cool Unix stuff
        cout << "LINUX detected\n";
        namespace fs = std::filesystem;
        if (!fs::is_directory("results") || !fs::exists("results")) { // Check if res folder exists
            fs::create_directory("results"); // create res folder
        }
    #else
        cout << "No supported os detected. To get results saved in a file create a folder named 'results'.\n";
    #endif 

    // clears output file
    ofstream outfile;
    outfile.open ("results/results.csv", ios::out | ios::trunc);
    outfile << 
    "n Particles;n Dimensions;n Metropolis Steps;Equilibration Fraction;Accepted Steps;Found Energy;Elapsed Time;n Parameters;Parameters (undefinedNumber)\n";
    outfile.close();
    outfile.open ("results/energies.csv", ios::out | ios::trunc);
    outfile.close();

    int methods[] = {1};

    time_point<system_clock> tot_time_start = high_resolution_clock::now();
    for (int nPar : numberOfParticles)
    {
        for (int nDim : numberOfDimensions)
        {
            for (double nAlpha : alpha)
            {
                for (int met : methods)
                {
                    System* system = new System(seed);
                    system->setHamiltonian              (new HarmonicOscillator(system, omega, true));
                    system->setWaveFunction             (new SimpleGaussian(system, nAlpha, dt));
                    system->setInitialState             (new RandomUniform(system, nDim, nPar));
                    system->setEquilibrationFraction    (equilibration);
                    system->setStepLength               (stepLength);
                    if (met == 0)
                    {
                        system->runMetropolisSteps          (numberOfSteps);
                        // cout << "Metropolis\n";
                    } else
                    {
                        system->runImportanceSamplingSteps  (numberOfSteps);
                        // cout << "Importance Sampling\n";
                    }
                }
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
