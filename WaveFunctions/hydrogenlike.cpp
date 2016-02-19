#include "hydrogenlike.h"
#include <cmath>
#include "system.h"
#include "particle.h"


HydrogenLike::HydrogenLike(System* system,
                           double alpha,
                           double beta) :
        WaveFunction(system) {
    m_numberOfParameters = 2;
    m_parameters = new double[2];
    m_parameters[0] = alpha;
    m_parameters[1] = beta;
}

double HydrogenLike::evaluate(Particle* particles) {
    double r1 = 0;
    double r2 = 0;
    double r12 = 0;

    for (int k=0; k<m_system->getNumberOfDimensions(); k++) {
        const double x1  = particles[0].getPosition()[k];
        const double x2  = particles[1].getPosition()[k];
        const double x12 = x2 - x1;
        r1  += x1  * x1;
        r2  += x2  * x2;
        r12 += x12 * x12;
    }
    r1  = sqrt(r1);
    r2  = sqrt(r2);
    r12 = sqrt(r12);

    const double interaction = m_parameters[1] < 1000 ? exp((0.5 * r12) / (1.0 + m_parameters[1] * r12)) : 1.0;
    const double single = exp(-m_parameters[0] * (r1 + r2));
    return single * interaction;
}

double HydrogenLike::computeKineticEnergy(Particle* particles) {
}
