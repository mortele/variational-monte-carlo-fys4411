#include "particle.h"
#include <cassert>

Particle::Particle(const std::vector<double> &position)
{
    m_numberOfDimensions = position.size();
    m_position = position;
    m_initialPosition = position; // Save initial position. Notice that this is a copy, not a reference.
    m_EquilibrationPosition = position;
}

void Particle::adjustPosition(double change, unsigned int dimension)
{
    m_position.at(dimension) += change;
}
void Particle::setPosition(double new_position, unsigned int dimension)
{
    m_position.at(dimension) = new_position;
}
void Particle::resetPosition()
{
    m_position = m_initialPosition;
}

void Particle::saveEquilibrationPosition()
{
    m_EquilibrationPosition = m_position;
}

void Particle::resetEquilibrationPosition()
{
    m_position = m_EquilibrationPosition;
}

double particle_r2(Particle &p)
{
    static const int numberOfDimensions = p.getNumberOfDimensions();
    double ret = 0;
    for (int q = 0; q < numberOfDimensions; q++)
    {
        ret += p.getPosition().at(q) * p.getPosition().at(q);
    }
    return ret;
}

double particle_r2(Particle &p1, Particle &p2)
{
    static const int numberOfDimensions = p1.getNumberOfDimensions();
    double rdiff;
    double ret = 0;
    for (int q = 0; q < numberOfDimensions; q++)
    {
        rdiff = p1.getPosition().at(q) - p2.getPosition().at(q);
        ret += rdiff * rdiff;
    }
    return ret;
}