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

#include <omp.h>

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


    int numberOfDimensions[]    = {1, 2, 3};
    // int numberOfParticles[]     = {1,10,100,500}; //{1, 2, 3}; 
    int numberOfParticles[]     = {1,10,100}; 
    int numberOfSteps           = (int) 1e4;
    double omega                = 1.0;              // Oscillator frequency.
    double alpha[]              = {.46}; // Variational parameter.
	// double alpha[] = {0.38, 0.42, 0.46, 0.5, 0.54, 0.58, 0.62};
    double stepLength           = 2;              // Metropolis step length.
    double equilibration        = 0.1;              // Amount of the total steps used
                                                    // for equilibration.
    int methods[]               = {0};
    double dt[]                 = {0.001};
    // double dt[]                 = {0.001, 0.005, 0.01};
    //for steepest descent
    bool do_steepest_descent    = true;
    double alpha_guess          = 0.45;
    int sd_steps                = (int) 1e3;
    int nIterations             = 1000;
    double eta                  = .001;
    double alphaChange          = 10;

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

    time_point<system_clock> tot_time_start = high_resolution_clock::now();

    // omp_set_num_threads(8);
    // cout << omp_get_num_threads() << endl;

    //checks if want to use steepest descent to optimize alpha
    if (do_steepest_descent)
    {    
        #pragma omp parallel for schedule(dynamic) default(shared) collapse(4)
        for (unsigned int nPar = 0; nPar < sizeof(numberOfParticles)/sizeof(numberOfParticles[0]); nPar++)
        {
            for (unsigned int nDim = 0; nDim < sizeof(numberOfDimensions)/sizeof(numberOfDimensions[0]); nDim++)
            {
                for (unsigned int met = 0; met < sizeof(methods)/sizeof(methods[0]); met++)
                {
                    for (unsigned int ddt = 0; ddt < sizeof(dt)/sizeof(dt[0]); ddt++)
                    {
                            
                        //steepest descent keeps goining until the desired number of iterations or until
                        //the change in alpha is acceptably small
                        int iters = nIterations;
                        for (int iter = 0; iter < nIterations; iter++)
                        {    
                            System* system = new System(seed);
                            system->setHamiltonian              (new HarmonicOscillator(system, omega, true));
                            system->setWaveFunction             (new SimpleGaussian(system, alpha_guess, dt[ddt]));
                            system->setInitialState             (new RandomUniform(system, numberOfDimensions[nDim], numberOfParticles[nPar]));
                            system->setEquilibrationFraction    (0); 
                            system->setStepLength               (stepLength);
                            
                            if (methods[met] == 0)
                            {
                                system->runMetropolisSteps         (sd_steps, false, "");
                            } else if (methods[met] == 1)
                            {
                                system->runImportanceSamplingSteps (sd_steps, false, "");
                            }
                            
                            double currEnergy = system->getSdRes()[0];
                            double currDeltaPsi = system->getSdRes()[1];
                            double currDerivativePsiE = system->getSdRes()[2];
                            alphaChange = eta*2*(currDerivativePsiE - currEnergy*currDeltaPsi);
                            alpha_guess -= alphaChange;

                            if (abs(alphaChange) < 1e-6 && abs(alpha_guess) < 2)
                            {
                                // cout << "iter: " << iter << endl;
                                iters = iter+1;
                                break;
                            }
                        }

                        string printStuff = "";
                        printStuff.append("Found best alpha: " + to_string(alpha_guess) + " after " + to_string(iters));
                        printStuff.append(" iterations on thread " + to_string(omp_get_thread_num()) + ".\n\n");


                        System* system = new System(seed);
                        system->setHamiltonian              (new HarmonicOscillator(system, omega, true));
                        system->setWaveFunction             (new SimpleGaussian(system, alpha_guess, dt[ddt]));
                        system->setInitialState             (new RandomUniform(system, numberOfDimensions[nDim], numberOfParticles[nPar]));
                        system->setEquilibrationFraction    (equilibration);
                        system->setStepLength               (stepLength);

                        if (methods[met] == 0)
                        {
                            system->runMetropolisSteps          (numberOfSteps, true, printStuff);
                        } else if (methods[met] == 1)
                        {
                            system->runImportanceSamplingSteps  (numberOfSteps, true, printStuff);
                        }
                    }
                }
            }
        }
    } else
    {
        #pragma omp parallel for schedule(dynamic) default(shared) collapse(5)
        for (unsigned int nPar = 0; nPar < sizeof(numberOfParticles)/sizeof(numberOfParticles[0]); nPar++)
        {
            for (unsigned int nDim = 0; nDim < sizeof(numberOfDimensions)/sizeof(numberOfDimensions[0]); nDim++)
            {
                for (unsigned int met = 0; met < sizeof(methods)/sizeof(methods[0]); met++)
                {
                    for (unsigned int nAlpha = 0; nAlpha < sizeof(alpha)/sizeof(alpha[0]); nAlpha++)
                    {
                        for (unsigned int ddt = 0; ddt < sizeof(dt)/sizeof(dt[0]); ddt++)
                        {
                            System* system = new System(seed);
                            system->setHamiltonian              (new HarmonicOscillator(system, omega, true));
                            system->setWaveFunction             (new SimpleGaussian(system, alpha[nAlpha], dt[ddt]));
                            system->setInitialState             (new RandomUniform(system, numberOfDimensions[nDim], numberOfParticles[nPar]));
                            system->setEquilibrationFraction    (equilibration);
                            system->setStepLength               (stepLength);
                            
                            if (methods[met] == 0)
                            {
                                system->runMetropolisSteps          (numberOfSteps, true, "");
                            } else if (methods[met] == 1)
                            {
                                system->runImportanceSamplingSteps  (numberOfSteps, true, "");
                            }
                        }
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
