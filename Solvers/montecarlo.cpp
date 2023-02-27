#include "montecarlo.h"
#include "Math/random.h"

MonteCarlo::MonteCarlo(std::unique_ptr<class Random> rng)
{
    m_rng = std::move(rng); // std::move transfers ownership of the object pointed to by rng to m_rng. Used to avoid copying large objects.
}
