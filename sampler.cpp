#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm> 
#include <sstream> 
#include <iterator> 
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

using namespace std;
using std::cout;
using std::endl;


Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
    m_numberOfSampledMetropolisSteps = steps - steps*m_system->getEquilibrationFraction();
}

void Sampler::sample(bool acceptedStep) {
    // Make sure the sampling variable(s) are initialized at the first step.
    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0;
        m_energies.clear();
        m_positions.clear();
    }

    /* Here you should sample all the interesting things you want to measure.
     * Note that there are (way) more than the single one here currently.
     */
    double localEnergy = m_system->getHamiltonian()->
                         computeLocalEnergy(m_system->getParticles());
    m_cumulativeEnergy  += localEnergy;
    m_energies.push_back(localEnergy);

    vector<class Particle*> particles = m_system->getParticles();
    vector<vector<double>> currPos = vector<vector<double>>();
    for (int i = 0; i < m_system->getNumberOfParticles(); i++)
    {
        currPos.push_back(particles[i]->getPosition());
    }
    m_positions.push_back(currPos);
    m_stepNumber++;
}

void Sampler::getOutput() {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibrationFraction();
    vector<double> pa = m_system->getWaveFunction()->getParameters();
    
    m_output.append("  -- System info -- \n");
    m_output.append(" Number of particles  : " + to_string(np) + "\n");
    m_output.append(" Number of dimensions : " + to_string(nd) + "\n");
    m_output.append(" Number of Metropolis steps run : 10^" + to_string(static_cast<int>(log10(ms))) + "\n");
    m_output.append(" Number of equilibration steps  : 10^" + to_string(static_cast<int>(log10(round(ms*ef)))) + "\n");
    // m_output.append(" Number of Metropolis steps run : 10^" + to_string(log10(ms)) + "\n");
    //m_output.append(" Number of equilibration steps  : 10^" + to_string(log10(round(ms*ef))) + "\n");
    m_output.append("\n");
    for (int i=0; i < p; i++) {
        m_output.append(" Parameter " + to_string(i+1) + " : " + to_string(pa.at(i)) + "\n");
    }
    m_output.append("\n");
    m_output.append("  -- Reults -- \n");
    m_output.append(" Energy : " + to_string(m_energy) + "\n");
    m_output.append("\n");
    
}

void Sampler::printOutputToTerminal() {
    
    cout << endl;
    cout << m_output;

    /* cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << log10(round(ms*ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
    }
    cout << endl;
    cout << "  -- Reults -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << endl; */
}

void Sampler::printOutputToFile() {
    
    ofstream outfile;
    outfile.open ("results.txt", ios::out | ios::trunc);
    outfile << m_output;
    outfile << endl;
    outfile << endl;

    /* ostringstream oss;
    copy(m_energies.begin(), m_energies.end()-1, ostream_iterator<vector<double>>(oss, ","));
    oss << m_energies.back();
    outfile << oss.str() << endl << endl; */

    /* copy(m_positions.begin(), m_positions.end()-1, ostream_iterator<double>(oss, ","));
    oss << m_energies.back();
    outfile << oss.str(); */
    outfile.close();

}

void Sampler::computeAverages() {
    /* Compute the averages of the sampled quantities. You need to think
     * thoroughly through what is written here currently; is this correct?
     */
    m_energy = m_cumulativeEnergy / (m_numberOfSampledMetropolisSteps);
}
