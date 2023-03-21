#include "particle.h"
#include <cassert>

Particle::Particle(const std::vector<double> &position)
{
    m_numberOfDimensions = position.size();
    m_position = position;
    m_initialPosition = position; // Save initial position. Notice that this is a copy, not a reference.
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

double particle_r2(Particle &p) 
{
    static const int numberOfDimensions = p.getNumberOfDimensions();
    double ret = 0;
    for(int q = 0; q < numberOfDimensions; q++) {
        ret += p.getPosition().at(q)*p.getPosition().at(q);
    }
    return ret;
}

double particle_r2(Particle &p1, Particle &p2)
{
    static const int numberOfDimensions = p1.getNumberOfDimensions();
    double rdiff;
    double ret = 0;
    for(int q = 0; q < numberOfDimensions; q++) {
        rdiff = p1.getPosition().at(q) - p2.getPosition().at(q);
        ret += rdiff*rdiff;
    }
    return ret;
}

double dot_product(std::vector<double> &v1, std::vector<double> &v2, int numberOfDimensions) 
{
    double ret = 0;
    for(int i = 0; i < numberOfDimensions; i++)
        ret += v1.at(i) * v2.at(i);

    return ret;
}

void particle_rdiff(std::vector<double> &diff, Particle &p1, Particle &p2, double scale = 1.0)
{
    static const int numberOfDimensions = p1.getNumberOfDimensions();
    for(int i = 0; i < numberOfDimensions; i++)
        diff.at(i) = scale*(p1.getPosition().at(i) - p2.getPosition().at(i));
}