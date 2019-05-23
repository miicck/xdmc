#ifndef __PARTICLE__
#define __PARTICLE__

// An abstract particle, which can interact with
// other particles (by default via the coulomb
// interaction).
class particle
{
public:
        static int count;
        particle();
        ~particle();

        virtual double interaction(particle* other); // The interaction energy with some other particle
        virtual void diffuse(double tau)=0;   // Called when a config diffuses in DMC
        virtual particle* copy()=0;           // Should return a (deep) copy of this particle
        virtual double charge()=0;            // The charge of this particle (electron charge = -1)
        virtual double mass()=0;              // The mass of this particle (electron mass = 1)
        virtual void sample_wavefunction()=0; // Called when a request is sent to sample a walker wvfn.

        // The location of this particle
        double* coords;
};

#endif
