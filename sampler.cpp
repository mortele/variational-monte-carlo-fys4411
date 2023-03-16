#include <memory>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "system.h"
#include "sampler.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

using std::fixed;
using std::setprecision;
using std::setw;

Sampler::Sampler(
    unsigned int numberOfParticles,
    unsigned int numberOfDimensions,
    double stepLength,
    unsigned int numberOfMetropolisSteps)
{
    m_stepNumber = 0;
    m_numberOfMetropolisSteps = numberOfMetropolisSteps;
    m_numberOfParticles = numberOfParticles;
    m_numberOfDimensions = numberOfDimensions;
    m_stepLength = stepLength;

    m_energy = 0;
    m_energy_variance = 0;
    m_energy_std = 0;

    m_cumulativeEnergy = 0;
    m_cumulativeEnergy2 = 0;
    m_numberOfAcceptedSteps = 0;

    m_numberOfParams = 1; // this should not be hard coded but we will change it later

    m_cumulativeDerPsiE = std::vector<double>(m_numberOfParams, 0);
    m_cumulativedeltaPsi = std::vector<double>(m_numberOfParams, 0);
    m_deltaPsi = std::vector<double>(m_numberOfParams, 0);
    m_derPsiE = std::vector<double>(m_numberOfParams, 0);
    m_energyDerivative = std::vector<double>(m_numberOfParams, 0);
}

void Sampler::sample(bool acceptedStep, System *system)
{
    /* Here you should sample all the interesting things you want to measure.
     * Note that there are (way) more than the single one here currently.
     */
    auto localEnergy = system->computeLocalEnergy();
    m_cumulativeEnergy += localEnergy;
    m_cumulativeEnergy2 += (localEnergy * localEnergy);
    m_stepNumber++;
    m_numberOfAcceptedSteps += acceptedStep;

    for (int i = 0; i < m_numberOfParams; i++)
    {
        m_deltaPsi.at(i) = system->computeParamDerivative(i);
        m_derPsiE.at(i) = m_deltaPsi.at(i) * localEnergy;
        m_cumulativedeltaPsi.at(i) += m_deltaPsi.at(i);
        m_cumulativeDerPsiE.at(i) += m_derPsiE.at(i);
    }
}

void Sampler::printOutputToTerminal(System &system)
{
    auto pa = system.getWaveFunctionParameters();
    auto p = pa.size();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << m_numberOfParticles << endl;
    cout << " Number of dimensions : " << m_numberOfDimensions << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(m_numberOfMetropolisSteps) << endl;
    cout << " Step length used : " << m_stepLength << endl;
    cout << " Ratio of accepted steps: " << ((double)m_numberOfAcceptedSteps) / ((double)m_numberOfMetropolisSteps) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (unsigned int i = 0; i < p; i++)
    {
        cout << " Parameter " << i + 1 << " : " << pa.at(i) << endl;
    }
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << " Energy variance : " << m_energy_variance << endl;
    cout << " Energy std : " << m_energy_std << endl;
    cout << endl;
}

void Sampler::writeOutToFile(System &system, std::string filename, double omega, bool analytical, bool importanceSampling)
{
    std::ifstream exsists_file(filename.c_str());

    std::fstream outfile;
    auto pa = system.getWaveFunctionParameters();
    int p = pa.size();
    int w = 20;

    if (!exsists_file.good())
    {
        outfile.open(filename, std::ios::out);
        outfile << setw(w) << "Dimensions"
                << setw(w) << "Particles"
                << setw(w) << "Metro-steps"
                << setw(w) << "Omega"
                << setw(w) << "StepLength";
        for (int i = 0; i < p; i++)
            outfile << setw(w - 1) << "WF" << (i + 1);

        outfile << setw(w) << "Energy"
                << setw(w) << "Energy_std"
                << setw(w) << "Energy_var"
                << setw(w) << "Accept_number"
                << setw(w) << "Accept_ratio"
                << setw(w) << "Imposampling"
                << setw(w) << "Analytical"
                << "\n";
    }
    else
    {
        outfile.open(filename, std::ios::out | std::ios::app);
    }
    outfile << setw(w) << m_numberOfDimensions
            << setw(w) << m_numberOfParticles
            << setw(w) << setprecision(5) << m_numberOfMetropolisSteps
            << setw(w) << fixed << setprecision(5) << omega
            << setw(w) << fixed << setprecision(5) << m_stepLength;

    for (int i = 0; i < p; i++)
        outfile << setw(w) << fixed << setprecision(5) << pa.at(i);

    outfile << setw(w) << fixed << setprecision(5) << m_energy
            << setw(w) << fixed << setprecision(5) << m_energy_std
            << setw(w) << fixed << setprecision(5) << m_energy_variance
            << setw(w) << fixed << setprecision(5) << m_numberOfAcceptedSteps
            << setw(w) << fixed << setprecision(5) << m_acceptRatio
            << setw(w) << fixed << setprecision(5) << importanceSampling
            << setw(w) << fixed << setprecision(5) << analytical
            << "\n";

    outfile.close();
}

void Sampler::output(System &system, std::string filename, double omega, bool analytical, bool importanceSampling)
{
    // Output information from the simulation, either as file or print
    if (filename == "") // this is dumbly duplicated now
    {
        printOutputToTerminal(system);
    }
    else
    {
        writeOutToFile(system, filename, omega, analytical, importanceSampling);
    }
}

void Sampler::computeAverages()
{
    /* Compute the averages of the sampled quantities.
     */
    m_energy = m_cumulativeEnergy / m_numberOfMetropolisSteps;

    m_cumulativeEnergy2 /= m_numberOfMetropolisSteps;
    m_energy_variance = (m_cumulativeEnergy2 - m_energy * m_energy);
    m_energy_std = sqrt(m_energy_variance);
    m_acceptRatio = ((double)m_numberOfAcceptedSteps) / ((double)m_numberOfMetropolisSteps);

    for (int i = 0; i < m_numberOfParams; i++)
    {
        m_derPsiE[i] = m_cumulativeDerPsiE[i] / m_numberOfMetropolisSteps;
        m_deltaPsi[i] = m_cumulativedeltaPsi[i] / m_numberOfMetropolisSteps;
    }
}

void Sampler::WriteTimingToFiles(System& system, std::string filename, bool analytical, unsigned int numberOfEquilibrationSteps, double timing){
    std::ifstream exsists_file(filename.c_str());

    std::fstream outfile;
    auto pa = system.getWaveFunctionParameters();
    int p = pa.size();
    int w = 20;

    if (!exsists_file.good())
    {
    
        outfile.open(filename, std::ios::out);
        outfile << setw(w) << "Dimensions"
                << setw(w) << "Particles"
                << setw(w) << "Metro-steps"
                << setw(w) << "Eq-steps"
                << setw(w) << "StepLength"
                << setw(w) << "Time"
                << setw(w) << "Analytical"
                << setw(w) << "Energy"
                << setw(w) << "Energy_std"
                << "\n";
    }
    else{
        outfile.open(filename, std::ios::out | std::ios::app);
    }
    outfile << setw(w) << m_numberOfDimensions
            << setw(w) << m_numberOfParticles
            << setw(w) << setprecision(5) << m_numberOfMetropolisSteps
            << setw(w) << setprecision(5) << numberOfEquilibrationSteps
            << setw(w) << fixed << setprecision(5) << m_stepLength
            << setw(w) << fixed << setprecision(0) << timing
            << setw(w) << analytical
            << setw(w) << setprecision(5) << m_energy
            << setw(w) << setprecision(5) << m_energy_std
            << "\n";



    outfile.close();

}


std::vector<double> Sampler::getEnergyDerivative()
{
    for (int i = 0; i < m_numberOfParams; i++)
    {
        m_energyDerivative[i] = 2 * (m_derPsiE[i] - m_energy * m_deltaPsi[i]);
    }
    return m_energyDerivative;
}
