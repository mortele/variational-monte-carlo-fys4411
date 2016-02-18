#pragma once
#include "harmonicoscillator.h"

class HarmonicOscillatorInteracting : public HarmonicOscillator {
public:
    HarmonicOscillatorInteracting(class System* system, double gamma);
    double computeLocalEnergy(Particle *particles);

private:
    double m_gamma  = 0;
    double m_a      = 0;
    double m_a2     = 0;
};
