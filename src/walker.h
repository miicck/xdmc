#ifndef __WALKER__
#define __WALKER__

#include <vector>

#include "particle.h"

// The object used by the diffusion monte carlo algorithm
// to represent a snapshot of the system.
class walker
{
public:
        static int count;
        double potential();
        void diffuse(double tau);
        walker();
        ~walker();
        walker* copy();
        void sample_wavefunction();

private:

        // The particles in this system snapshot
        std::vector<particle*> particles;

        // Create a walker from a given particle set 
        walker(std::vector<particle*> particles);
};

#endif
