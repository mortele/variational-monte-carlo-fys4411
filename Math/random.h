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

