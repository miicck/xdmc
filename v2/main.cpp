#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>

// Constants
const double PI = 3.141592653589793;

// Program parameters
int target_population = 1000;
int dmc_itterations   = 10000;
double tau            = 0.01;
double trial_energy   = 0;

// Generate a uniform random number \in [0,1]
double rand_uniform()
{
	return (rand()%RAND_MAX)/double(RAND_MAX);
}

// Generate a normal random number with variance var
// using a box-muller transformation
double rand_normal(double var)
{
	double u1 = rand_uniform();
	double u2 = rand_uniform();
	return sqrt(-2*var*log(u1)) * sin(2*PI*u2);
}

// An abstract particle
class particle
{
public:
	virtual double interaction(particle* other)
	{
		double dx = this->x - other->x;
		double dy = this->y - other->y;
		double dz = this->z - other->z;
		double r = sqrt(dx*dx+dy*dy+dz*dz);
		return (this->charge()*other->charge())/r;
	}

	virtual void diffuse(double tau)=0;
	virtual particle* copy()=0;
	virtual double charge()=0;
	virtual double mass()=0;
	double x;
	double y;
	double z;
};

// A particle that is described as point-like
// and static within the DMC algorithm
class classical_particle : public particle
{
public:
	virtual void diffuse(double tau) { }
};

// A particle who is described by an ensemble of
// diffusing walkers within the DMC algorithm
class quantum_particle : public particle
{
public:
	virtual void diffuse(double tau)
	{
		x += rand_normal(tau);
		y += rand_normal(tau);
		z += rand_normal(tau);
	}
};

// e-
class electron : public quantum_particle
{
	virtual double charge() { return -1; }
	virtual double mass()   { return  1; }

	virtual particle* copy()
	{
		electron* ret = new electron();
		ret->x = x;
		ret->y = y;
		ret->z = z;
		return ret;
	}
};

// An atomic nucleus
class nucleus : public classical_particle
{
public:
	nucleus(double charge, double mass)
	{
		atomic_number = charge;
		mass_number   = mass;
	}
	virtual double charge() { return atomic_number; }
	virtual double mass()   { return mass_number;   }

	virtual particle* copy()
	{
		nucleus* ret = new nucleus(atomic_number, mass_number);
		ret->x = x;
		ret->y = y;
		ret->z = z;
		return ret;
	}
private:
	double atomic_number;
	double mass_number;
};

// A snapshot of a system used to describe
// it via diffusion monte carlo
class walker
{
public:
	walker(std::vector<particle*> particles)
	{	
		this->particles = particles;
	}

	~walker()
	{
		for (int i=0; i<particles.size(); ++i)
			delete(particles[i]);
	}

	double potential()
	{
		double ret = 0;
		for (int i = 0; i < particles.size(); ++i)
		{
			particle* p1 = particles[i];
			for (int j=0; j<i; ++j)
			{
				particle* p2 = particles[j];
				ret += p1->interaction(p2);
			}
		}
		return ret;
	}

	void diffuse(double tau)
	{
		for (int i=0; i<particles.size(); ++i)
		{
			particle* p = particles[i];
			p->diffuse(tau);
		}
	}

	walker* copy()
	{
		std::vector<particle*> copied_particles;
		for (int i=0; i<particles.size(); ++i)
			copied_particles.push_back(particles[i]->copy());
		return new walker(copied_particles);
	}

private:
	std::vector<particle*> particles;
};

// Returns a set of particles describing our system
std::vector<particle*> get_system()
{
	std::vector<particle*> ret;
	ret.push_back(new electron());
	nucleus* n = new nucleus(1,1);
	n->x = 0;
	n->y = 0;
	n->z = 0;
	ret.push_back(n);
	return ret;
}

// Returns the number of walkers that should survive after a 
// walker moves from potential pot_before to pot_after. i.e returns
// 0 if the walker should die and 1 if it should live, 2 if another should
// spawn etc...
int walkers_surviving(double pot_before, double pot_after)
{
	double p = exp(-tau*(pot_before + pot_after - 2*trial_energy)/2);
	int ret = 0;
	while(p>1)
	{
		++ret;
		--p;
	}
	if (rand_uniform() < p)
		++ret;
	return ret;
}

// Contains the results of a DMC iteration
struct iter_result
{
	int population;
	double average_potential;

	void print()
	{
		std::cout << population << ", " << average_potential << "\n";
	}

	static std::ofstream file;
	void write()
	{
		file << population << "," << average_potential << "\n";
	}
};
std::ofstream iter_result::file("out");

// The main DMC algorithm
int main(int argc, char** argv)
{
	std::vector<walker*> walkers;
	for (int i=0; i<target_population; ++i)
	{
		walker* w = new walker(get_system());
		w->diffuse(1.0);
		walkers.push_back(w);
	}
	
	for (int iter = 1; iter<dmc_itterations+1; ++iter)
	{
		iter_result ir;

		std::vector<walker*> new_walkers;
		for (int n=0; n<walkers.size(); ++n)
		{
			walker* w = walkers[n];
			double pot_before = w->potential();
			w->diffuse(tau);
			double pot_after  = w->potential();
			int surviving = walkers_surviving(pot_before, pot_after);

			ir.average_potential += surviving*(pot_before+pot_after)/2;

			for (int s=0; s<surviving; ++s)
				new_walkers.push_back(w->copy());
			delete(w);
		}

		walkers = new_walkers;

		ir.population = walkers.size();
		ir.average_potential /= ir.population;

		trial_energy = ir.average_potential;
		
		ir.write();
	}
}
