#ifndef __PARTICLE__
#define __PARTICLE__

#include "constants.h"

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

// A particle that is described as point-like
// and static within the DMC algorithm.
class classical_particle : public particle
{
public:
        // Classical particles don't diffuse
        virtual void diffuse(double tau) { }
        virtual void sample_wavefunction() { }
};

// A particle who is described by an ensemble of
// diffusing walkers within the DMC algorithm
class quantum_particle : public particle
{
public:
        virtual void diffuse(double tau);
        virtual void sample_wavefunction();
};

// I read about these in a physics textbook once, thought
// they might be important.
class electron : public quantum_particle
{
        virtual double charge() { return -1; }
        virtual double mass()   { return  1; }
        virtual particle* copy();
};

// A fermion with no interactions
class non_interacting_fermion : public quantum_particle
{
public:
        virtual double interaction(particle* other) { return 0; }
        virtual double mass()   { return 1; }
        virtual double charge() { return 0; }
        particle* copy();
};

// An atomic nucleus, as described by an atomic
// number and a mass number.
class nucleus : public classical_particle
{
public:
        nucleus(double atomic_number, double mass_number);
        virtual double charge() { return atomic_number; }
        virtual double mass()   { return mass_number * AMU; } // Convert mass to atomic units
        virtual particle* copy();
private:
        double atomic_number;
        double mass_number;
};

#endif
