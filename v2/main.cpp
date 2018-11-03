#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <mpi.h>

// Constants
const double PI = 3.141592653589793;

// Program parameters
int target_population = 500;
int dmc_iterations    = 10000;
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

// An abstract particle, which can interact with
// other particles (by default using the coulomb
// interaction)
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
		// Diffuse the particle by moving each
		// coordinate by an amount sampled from
		// a normal distribution with variance tau
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
		// Copy this electron
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
		// Copy this nucleus
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

// The object used by the diffusion monte carlo algorithm
// to represent a snapshot of the system.
class walker
{
public:
	walker(std::vector<particle*> particles)
	{	
		// Setup the walker with the particles
		// that describe the system.
		this->particles = particles;
	}

	~walker()
	{
		// Clear up memory (delete all the particles).
		for (int i=0; i<particles.size(); ++i)
			delete(particles[i]);
	}

	double potential()
	{
		// Evaluate the potential of the system
		// in the configuration described by this
		// walker; this is a sum of particle-particle
		// interactions, avoiding double counting.
		double ret = 0;
		for (int i = 0; i < particles.size(); ++i)
			for (int j=0; j<i; ++j) // Note j<i => no double counting
				ret += particles[i]->interaction(particles[j]);
		return ret;
	}

	void diffuse(double tau)
	{
		// Diffuse all of the particles (classical
		// particles will automatically not diffuse)
		for (int i=0; i<particles.size(); ++i)
			particles[i]->diffuse(tau);
	}

	walker* copy()
	{
		// Return a copy of this walker (copy each of the particles)
		std::vector<particle*> copied_particles;
		for (int i=0; i<particles.size(); ++i)
			copied_particles.push_back(particles[i]->copy());
		return new walker(copied_particles);
	}

private:
	
	// The particles in this system snapshot
	std::vector<particle*> particles;
};

// Returns the set of particles in our system
std::vector<particle*> get_system()
{
	std::vector<particle*> ret;

	// For now, a hard coded hydrogen atom
	ret.push_back(new electron());
	nucleus* n = new nucleus(1,1);
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
	int ret = int(floor(p));
	p -= ret;
	if (rand_uniform() < p)
		++ret;
	return ret;
}

// Contains the results of a single DMC iteration
struct iter_result
{
	int population;
	double average_potential;

	// Directly print the results 
	void print()
	{
		std::cout << population << ", " << average_potential << "\n";
	}
	
	// Write the resuls to the "out" file
	void write()
	{
		static std::ofstream file("out");
		file << population << "," << average_potential << "\n";
	}

	// Reduce the results of an iteration across processes
	void mpi_reduce(int np, int pid)
	{
		double av_pot = 0;
		int pop_sum = 0;
		int ierr    = 0;

		ierr = MPI_Reduce(&population, &pop_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		ierr = MPI_Reduce(&average_potential, &av_pot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		if (pid == 0)
		{
			population = pop_sum;
			average_potential = av_pot / double(np);
		}
	}
};

// The main DMC algorithm
int main(int argc, char** argv)
{
	// Initialize mpi
	int ierr = MPI_Init(&argc, &argv);
	if (ierr != 0)
	{
		std::cout << "Fatal error, could not initialize MPI.";
		exit(1);
	}

	// Get the number of processes and my id within them
	int np = 1;
	int pid = 0;
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &np);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	srand(pid*clock());
	std::cout << "Process " << (pid+1) << "/" << np << " reporting for duty.\n";

	// Our DMC walkers
	std::vector<walker*> walkers;

	// Initialize the walkers
	for (int i=0; i<target_population; ++i)
	{
		walker* w = new walker(get_system());
		walkers.push_back(w);

		// Carry out an initial large diffusion
		// to avoid unneccasary equilibriation
		// and exact particle overlap at the origin
		w->diffuse(1.0);
	}
	
	// Run our DMC iterations
	for (int iter = 1; iter <= dmc_iterations; ++iter)
	{
		iter_result ir;

		// The new walkers that appear after the iteration
		std::vector<walker*> new_walkers;
		for (int n=0; n<walkers.size(); ++n)
		{
			walker* w = walkers[n];

			// Diffuse the walker and evaluate potential change
			double pot_before = w->potential();
			w->diffuse(tau);
			double pot_after  = w->potential();

			// Work out how many survivded the move
			int surviving = walkers_surviving(pot_before, pot_after);

			if (surviving == 0) delete(w); // Delete this walker if it died
			else new_walkers.push_back(w); // Otherwise it survives

			// Spawn new walkers
			for (int s=0; s < surviving-1; ++s) 
				new_walkers.push_back(w->copy());

			// Accumulate the average potential
			ir.average_potential += surviving*(pot_before+pot_after)/2;
		}

		// Set our walkers to the newly generated ones
		walkers = new_walkers;

		// Record the resuts of the iteration
		ir.population = walkers.size();
		ir.average_potential /= ir.population;
		ir.mpi_reduce(np, pid);
		if (pid == 0) ir.write();

		// E_T = <potential> to avoid population explosion/collapse
		trial_energy = ir.average_potential; 
	}

	// End mpi safely
	MPI_Finalize();
}
