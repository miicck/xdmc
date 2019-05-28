/*
    DMCPP
    Copyright (C) 2019 Michael Hutcheon (email mjh261@cam.ac.uk)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    For a copy of the GNU General Public License see <https://www.gnu.org/licenses/>.
*/

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

   	// -1 for identical fermions, 1 for identical bosons, 0 otherwise
	virtual int exchange_symmetry(particle* other)=0;

        // The location of this particle
        double* coords;

	// The probability that this particle
	// will diffuse into the given copy of it
	double overlap_prob(particle* clone);
};

// A particle that is described as point-like
// and static within the DMC algorithm.
class classical_particle : public particle
{
public:
        // Classical particles don't diffuse
        virtual void diffuse(double tau) { }
        virtual void sample_wavefunction() { }
	virtual int exchange_symmetry(particle* other) { return 0; }
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

	virtual int exchange_symmetry(particle* other)
	{ 
		if (dynamic_cast<electron*>(other))
			return -1; 
		return 0;
	}
};

// A fermion with no interactions
class non_interacting_fermion : public quantum_particle
{
public:
        virtual double interaction(particle* other) { return 0; }
        virtual double mass()   { return 1; }
        virtual double charge() { return 0; }
        particle* copy();

	virtual int exchange_symmetry(particle* other)
	{
		if (dynamic_cast<non_interacting_fermion*>(other))
			return -1;
		return 0;
	}
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


