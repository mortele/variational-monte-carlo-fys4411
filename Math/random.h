/*
 This is based on the ran1-generator from lib.cpp (see http://www.uio.no/studier/emner/matnat/fys/FYS3150/h14/index.html )
     ** The function
     **           ran1()
     ** is an "Minimal" random number generator of Park and Miller
     ** (see Numerical recipe page 280) with Bays-Durham shuffle and
     ** added safeguards. Call with idum a negative integer to initialize;
     ** thereafter, do not alter idum between sucessive deviates in a
     ** sequence. RNMX should approximate the largest floating point value
     ** that is less than 1.
     ** The function returns a uniform deviate between 0.0 and 1.0
     ** (exclusive of end-point values).
*/

#pragma once

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

class Random {
public:
    static long iy;
    static long iv[NTAB];
    static long seed;
    static int    nextInt(int upperLimit);
    static double nextDouble();
    static double nextGaussian(double mean, double standardDeviation);
    static void setSeed(long seed);
};

