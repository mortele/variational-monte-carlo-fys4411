#include "montecarlo.h"
#include "Math/random.h"

MonteCarlo::MonteCarlo(std::unique_ptr<class Random> rng)
{
    m_rng = std::move(rng);
}
