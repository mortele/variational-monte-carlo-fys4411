#include "random.h"
#include <cmath>

long     Random::iy = 0;
long     Random::iv[NTAB];
long     Random::seed = -1;
void Random::setSeed(long seed) {
    Random::seed = seed;
}

double Random::nextGaussian(double mean, double standardDeviation) {
    double standardNormalRandomNumber = sqrt( -2.0*log(1.0 - nextDouble()) ) * cos( 6.283185307 * nextDouble() );
    return standardDeviation*standardNormalRandomNumber + mean;
}

int Random::nextInt(int upperLimit) {
    return std::floor(nextDouble() * upperLimit);
}

double Random::nextDouble()
{
    int             j;
    long            k;
    double          temp;
    if (Random::seed <= 0 || !iy) {
        if (-(Random::seed) < 1) Random::seed=1;
        else Random::seed = -(Random::seed);
        for(j = NTAB + 7; j >= 0; j--) {
            k     = (Random::seed)/IQ;
            Random::seed = IA*(Random::seed - k*IQ) - IR*k;
            if(Random::seed < 0) Random::seed += IM;
            if(j < NTAB) iv[j] = Random::seed;
        }
        iy = iv[0];
    }
    k     = (Random::seed)/IQ;
    Random::seed = IA*(Random::seed - k*IQ) - IR*k;
    if(Random::seed < 0) Random::seed += IM;
    j     = iy/NDIV;
    iy    = iv[j];
    iv[j] = Random::seed;
    if((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

