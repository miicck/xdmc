#ifndef __POTENTIAL__
#define __POTENTIAL__

#include "particle.h"

class external_potential
{
public:
        virtual double potential(particle* p)=0;
};

#endif
