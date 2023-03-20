 #include <memory>
#include <cmath>
#include <cassert>

#include "interactinggaussian.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>

InteractingGaussian::InteractingGaussian(double alpha, double a, int num_particles) // a is the parameter in the interaction
{
    assert(alpha >= 0);
    m_interactionTerm = a;               // arameter in the interaction
    m_numberOfParameters = 1;            // this should not be hard coded
    m_numberOfParticles = num_particles; // this is the number of particles in the system. Different from the non interacting case, it is important
                                         // that this is a member variable, since we need to know the number of particles in the system in order to
                                         // compute the interaction term in the quantum force, for example, which do not currently have access to.
    m_parameters.reserve(1);
    m_parameters.push_back(alpha); // m_parameters is the vector of variational parameters
}

double InteractingGaussian::evaluate(std::vector<std::unique_ptr<class Particle>> &particles)
{
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i]->getPosition()
     * function.
     */
    int numberOfDimensions = particles.at(0)->getNumberOfDimensions();

    double r2 = 0;  // r2 is the sum of the squared coordinates of the r vector
    double r_q = 0; // r_q is the q'th coordinate of the r vector
    double alpha = m_parameters.at(0);
    double a = m_interactionTerm; // renaming for simplicity in formulas


    for (int i = 0; i < m_numberOfParticles; i++)
    {
        Particle particle = *particles.at(i);
        for (int q = 0; q < numberOfDimensions; q++)
        {
            r_q = particle.getPosition().at(q);
            r2 += r_q * r_q;
        }
    }

    double gaussian = std::exp(-alpha * r2); // this is the gaussian part of the wave function
    
    // now we have interaction, so instead of just the gaussian, we have to multiply by the interaction term
    double interaction = 1;
    double r_ij_q = 0;   // r_ij_q is the q'th coordinate of the r_ij vector
    double r_ij = 0;     // r_ij is the squared distance between particles i and j. This is important now that we have interaction
    double norm_rij = 0; // norm_r_ij is the distance between particles i and j
    
    
    for (int i = 0; i < m_numberOfParticles; i++)
    {
        Particle particle_i = *particles.at(i);
        for (int j = i + 1; j < m_numberOfParticles; j++) // we only need to loop over the particles with higher index than i because of the symmetry of the interaction
        {
            r_ij = 0;
            Particle particle_j = *particles.at(j);
            for (int q = 0; q < numberOfDimensions; q++)
            {
                r_ij_q = particle_i.getPosition().at(q) - particle_j.getPosition().at(q); // ex: r_ij_x = r_i_x - r_j_x
                r_ij += r_ij_q * r_ij_q;                                                  // ex: r_ij = r_ij_x^2 + r_ij_y^2 + r_ij_z^2
            }
            norm_rij = std::sqrt(r_ij); // there might be a more efficient way to do this, using u = ln(f(r_ij)).
            if (norm_rij > a)
            {
                interaction *= (1.0 - a / norm_rij);
            }
        }
    }

    return gaussian * interaction;
}

double InteractingGaussian::computeParamDerivative(std::vector<std::unique_ptr<class Particle>> &particles, int parameterIndex)
{
    /* Note that by derivative, we actually
     * mean the derivative with respect to the variational parameters.
     * Thanks to the good God, the interaction term does not depend on the variational parameters,
     * so we can just return the derivative of the gaussian part wrt alpha.
     * Notice we don't even need to multiply by the interaction term, because we then divide by the WF (TRIPLE CHECK THIS!).
     */
    int numberOfDimensions = particles.at(0)->getNumberOfDimensions();
    double parameter = m_parameters.at(parameterIndex); // this is not used now, but can be used when we generalize

    double r2_sum = 0;
    double r_q = 0;

    for (int k = 0; k < m_numberOfParticles; k++)
    {
        Particle particle = *particles.at(k);
        for (int q = 0; q < numberOfDimensions; q++)
        {
            r_q = particle.getPosition().at(q);
            r2_sum += r_q * r_q;
        }
    }
    return -r2_sum; // analytic derivative wrt alpha, only 1 param to optimize now, This needs to be generalized
}

// I still need to do this by hand to figure out the terms given our wf.
double InteractingGaussian::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>> &particles)
{
    /*All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle. And also we devide by the wave function because it simplifies the
     * calculations in the Metropolis algorithm and it is all we actually need.
     */
    int numberOfDimensions = particles.at(0)->getNumberOfDimensions();
    double alpha = m_parameters.at(0);
    double a = m_interactionTerm; // renaming for simplicity in formulas

    double r2_sum = 0;
    double r_q = 0;

    // first we compute the double derivative of the gaussian part
    for (int k = 0; k < m_numberOfParticles; k++)
    {
        Particle particle = *particles.at(k);
        for (int q = 0; q < numberOfDimensions; q++)
        {
            r_q = particle.getPosition().at(q);
            r2_sum += r_q * r_q;
        }
    }

    double gaussian_double_derivative = 2 * alpha * (2 * alpha * r2_sum - m_numberOfParticles * numberOfDimensions); // analytic double derivative

    // now we compute the double derivative of the interaction part
    double interaction_double_derivative = 0; // this will contain 3 more terms * ONLY if r_ij > a*

    double term1 = 0;
    double term2 = 0;
    double term3 = 0;

    for (int k = 0; k < m_numberOfParticles; k++) // we can join the loops later
    {
        Particle particle_k = *particles.at(k);
        // loop over all particles but not particle k
        for (int i = 0; i < m_numberOfParticles; i++)
        {

            Particle particle_i = *particles.at(i);
            double r_ik = 0;
            double r_ik_q = 0;
            double r_ik_norm = 0;
            double r_jk_norm = 0;
            double r_jk_r_ik = 0;
            double r_ik_q_r_k_q = 0; //  q'th coordinate of the r_ik * r_k vector
            double r_ik_rk = 0;      // r_ik * r_k
            for (int q = 0; q < numberOfDimensions; q++)
            {
                r_ik_q = particle_k.getPosition().at(q) - particle_i.getPosition().at(q);
                r_ik += r_ik_q * r_ik_q;

                r_ik_q_r_k_q = r_ik_q * r_ik_q * particle_k.getPosition().at(q);
                r_ik_rk += r_ik_q_r_k_q * r_ik_q_r_k_q;
            }
            r_ik_norm = std::sqrt(r_ik);
            if (r_ik_norm > a) // notice this eliminates the self-interaction term
            {
                term1 += r_ik_rk / (r_ik_norm * r_ik_norm * (r_ik_norm - a));
                term2 += 1 / (r_ik_norm * r_ik_norm) * (r_ik_norm - a) * (r_ik_norm - a);

                for (int j = 0; j < m_numberOfParticles; j++)
                {
                    Particle particle_j = *particles.at(j);
                    double r_jk = 0;
                    double r_jk_q = 0;
                    for (int q = 0; q < numberOfDimensions; q++)
                    {
                        r_jk_q = particle_k.getPosition().at(q) - particle_j.getPosition().at(q);
                        r_jk += r_jk_q * r_jk_q;

                        r_jk_r_ik += r_jk_q * r_ik_q; // this is the dot product of r_jk and r_ik
                    }
                    r_jk_norm = std::sqrt(r_jk);
                    if (r_jk_norm > a)
                    {
                        term3 += r_jk_r_ik / (r_ik_norm * r_jk_norm * (r_ik_norm - a) * (r_jk_norm - a));
                    }
                }
            }
        }
    }
    interaction_double_derivative = -4 * alpha * a * term1 - a * a * term2 + a * a * term3;

    return gaussian_double_derivative + interaction_double_derivative;
}

double InteractingGaussian::evaluate_w(int proposed_particle_idx, class Particle &proposed_particle, class Particle &old_particle, std::vector<std::unique_ptr<class Particle>> &particles)
{
    /*
     This is the wave function ratio for the Metropolis algorithm.
     It is a clever way to avoid having to evaluate the wave function for all particles at each step.
     Notice that the interaction term does not depend on the proposed particle, so we can just evaluate the gaussian part (I THINK).
    */
    static const int numberOfDimensions = particles.at(0)->getNumberOfDimensions(); // static to avoid redeclaration between calls
    static const double a = m_interactionTerm;
    const double alpha = m_parameters.at(0);

    double r2_proposed, r2_old;
    r2_proposed = 0;
    r2_old = 0;

    r2_proposed = particle_r2(proposed_particle);
    r2_old = particle_r2(old_particle);

    double gaussian = std::exp(-2.0 * alpha * (r2_proposed - r2_old));

    double interaction = 1;
    double r_gj_prime = 0;
    double r_gj = 0;
    double delta = 0;
    for(int i = 0; i < proposed_particle_idx; i++) 
    {
        r_gj_prime = std::sqrt( particle_r2(proposed_particle, *particles.at(i)) ); 
        r_gj = std::sqrt( particle_r2(old_particle, *particles.at(i)) );
        delta = (r_gj_prime > a)*(r_gj > a);
        interaction *= (1.0 - a/r_gj_prime) / (1.0 - a/r_gj);
    }
    for(int i = proposed_particle_idx+1; i < m_numberOfParameters; i++) {
        r_gj_prime = std::sqrt( particle_r2(proposed_particle, *particles.at(i)) ); 
        r_gj = std::sqrt( particle_r2(old_particle, *particles.at(i)) );
        delta = (r_gj_prime > a)*(r_gj > a);
        interaction *= (1.0- a/r_gj_prime) / (1.0- a/r_gj);
    }

    return gaussian * interaction * interaction;
}

// Notice that now we need to pass the whole vector of particles, because we need to compute the interaction term.
void InteractingGaussian::quantumForce(std::vector<std::unique_ptr<class Particle>> &particles, Particle &particle, std::vector<double> &force)
{
    static const int numberOfDimensions = particle.getNumberOfDimensions(); // static to avoid redeclaration between calls
    static const double a = m_interactionTerm;
    const double alpha = m_parameters.at(0);

    for (int q = 0; q < numberOfDimensions; q++)
    {
        force.at(q) = -4.0 * alpha * particle.getPosition().at(q); // analytic derivative wrt r_q for the gaussian part
    }

    // Now we need to add the interaction term
    double r_ij, r_ij_q, norm_rij;
    for (int i = 0; i < m_numberOfParticles; i++)
    {
        Particle &other_particle = *particles.at(i);
        norm_rij = 0;
        for (int q = 0; q < numberOfDimensions; q++)
        {
            r_ij_q = particle.getPosition().at(q) - other_particle.getPosition().at(q);
            r_ij += r_ij_q * r_ij_q;
        }
        norm_rij = std::sqrt(r_ij);
        if (norm_rij > a) // notice this automatically excludes the particle itself
        {
            for (int q = 0; q < numberOfDimensions; q++)
            {
                r_ij_q = particle.getPosition().at(q) - other_particle.getPosition().at(q); // NOTICE THIS IS INEFFICIENT, WE ALREADY COMPUTED IT BUT WAS NOT STORED
                force.at(q) *= a * r_ij_q / (norm_rij * norm_rij * (norm_rij - a));
            }
        }
    }
}

void InteractingGaussian::setParameters(std::vector<double> parameters)
{
    assert((int)parameters.size() == m_numberOfParameters);
    m_parameters = parameters;
}