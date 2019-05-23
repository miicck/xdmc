#ifndef __POTENTIAL__
#define __POTENTIAL__

#include "particle.h"

class external_potential
{
public:
        virtual double potential(particle* p)=0;
};

class harmonic_well : public external_potential
{
public:
        harmonic_well(double omega) { this->omega = omega; }
        virtual double potential(particle* p);
private:
        double omega = 1;
};

#endif
