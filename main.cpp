#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iterator>
#include <mpi.h>
#include <unistd.h>

// Constants
const double PI = 3.141592653589793; // You should know this one
const double AMU = 1822.88486193;    // Unified atomic mass unit / electron mass
const int MPI_ERROR = 1;	     // Thrown when an MPI call fails

// Forward declerations
class particle;
class external_potential;

// Program parameters
int dimensions        = 3;      // The dimensionality of our system
int target_population = 500;	// The number of DMC walkers per process
int dmc_iterations    = 1000;   // The number of DMC iterations to carry out
double tau            = 0.01;   // The DMC timestep
double trial_energy   = 0;      // Energy used to control the DMC population

// The system which will be copied to generate walkers
std::vector<particle*> template_system;
std::vector<external_potential*> potentials;

// Output files
std::ofstream wavefunction_file;
std::ofstream out_file;

// Generate a uniform random number \in [0,1].
double rand_uniform()
{
	return (rand()%RAND_MAX)/double(RAND_MAX);
}

// Generate a zero-mean noramlly distributed number
// with the specified variance
double rand_normal(double var)
{
	double u1 = rand_uniform();
	double u2 = rand_uniform();
	return sqrt(-2*var*log(u1)) * sin(2*PI*u2);
}

// An abstract particle, which can interact with
// other particles (by default via the coulomb
// interaction).
class particle
{
public:
	static int count;
	particle()
	{
		++ count;
		coords = new double[dimensions];
	}

	~particle()
	{
		-- count;
		delete[] coords;
	}

	virtual double interaction(particle* other)
	{
		// Coulomb
		double r = 0;
		for (int i=0; i<dimensions; ++i)
		{
			double dxi = this->coords[i] - other->coords[i];
			r += dxi * dxi;
		}
		r = sqrt(r);
		return (this->charge()*other->charge())/r;
	}

	virtual void diffuse(double tau)=0;   // Called when a config diffuses in DMC
	virtual particle* copy()=0;	      // Should return a (deep) copy of this particle
	virtual double charge()=0;	      // The charge of this particle (electron charge = -1)
	virtual double mass()=0;	      // The mass of this particle (electron mass = 1)
	virtual void sample_wavefunction()=0; // Called when a request is sent to sample a walker wvfn.

	// The location of this particle
	double* coords;
};
int particle::count = 0;

class external_potential
{
public:
	virtual double potential(particle* p)=0;
};

class harmonic_well : public external_potential
{
public:
	harmonic_well(double omega) { this->omega = omega; }

	virtual double potential(particle* p)
	{
		double r2 = 0;
		for (int i=0; i<dimensions; ++i)
			r2 += p->coords[i]*p->coords[i];
		return 0.5*r2*omega*omega;
	}

private:
	double omega = 1;
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
	virtual void diffuse(double tau)
	{
		// Diffuse the particle by moving each
		// coordinate by an amount sampled from
		// a normal distribution with variance tau.
		for (int i=0; i<dimensions; ++i)
			this->coords[i] += rand_normal(tau);
	}

	virtual void sample_wavefunction()
	{
		for (int  i=0; i<dimensions-1; ++i)
			wavefunction_file << this->coords[i] << ",";
		wavefunction_file << this->coords[dimensions-1];
	}
};

// I read about these in a physics textbook once, thought
// they might be important.
class electron : public quantum_particle
{
	virtual double charge() { return -1; }
	virtual double mass()   { return  1; }

	virtual particle* copy()
	{
		// Copy this electron.
		electron* ret = new electron();
		for (int i=0; i<dimensions; ++i)
			ret->coords[i] = this->coords[i];
		return ret;
	}
};

// A fermion with no interactions
class non_interacting_fermion : public quantum_particle
{
public:
	virtual double interaction(particle* other) { return 0; }
	virtual double mass()   { return 1; }
	virtual double charge() { return 0; }

	particle* copy()
	{
		particle* ret = new non_interacting_fermion();
		for (int i=0; i<dimensions; ++i)
			ret->coords[i] = this->coords[i];
		return ret;
	}
};

// An atomic nucleus, as described by an atomic
// number and a mass number. 
class nucleus : public classical_particle
{
public:
	nucleus(double atomic_number, double mass_number)
	{
		this->atomic_number = atomic_number;
		this->mass_number   = mass_number;
	}
	virtual double charge() { return atomic_number; }
	virtual double mass()   { return mass_number * AMU; } // Convert mass to atomic units

	virtual particle* copy()
	{
		// Copy this nucleus
		nucleus* ret = new nucleus(atomic_number, mass_number);
		for (int i=0; i<dimensions; ++i)
			ret->coords[i] = this->coords[i];
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
	static int count;

	double potential()
	{
		// Evaluate the potential of the system
		// in the configuration described by this walker
		double ret = 0;
		for (int i = 0; i < particles.size(); ++i)
		{
			// Sum up external potential contributions
			for (int j=0; j<potentials.size(); ++j)
				ret += potentials[j]->potential(particles[i]);

			// Particle-particle interactions
			// note j<i => no double counting
			for (int j=0; j<i; ++j)
				ret += particles[i]->interaction(particles[j]);
		}
		return ret;
	}

	void diffuse(double tau)
	{
		// Diffuse all of the particles (classical
		// particles will automatically not diffuse)
		for (int i=0; i<particles.size(); ++i)
			particles[i]->diffuse(tau);
	}

	walker()
	{	
		// Setup the walker with the particles
		// that describe the system.
		++ count;
		this->particles = create_system_particles();
	}

	~walker()
	{
		// Clear up memory (delete all the particles).
		-- count;
		for (int i=0; i<particles.size(); ++i)
			delete(particles[i]);
	}

	walker* copy()
	{
		// Return a copy of this walker (copy each of the particles)
		std::vector<particle*> copied_particles;
		for (int i=0; i<particles.size(); ++i)
			copied_particles.push_back(particles[i]->copy());
		return new walker(copied_particles);
	}

	void sample_wavefunction()
	{
		for (int i=0; i<particles.size(); ++i)
		{
			particles[i]->sample_wavefunction();
			if (i < particles.size() - 1)
				wavefunction_file << ";";
		}
		wavefunction_file << "\n";
	}

private:
	
	// The particles in this system snapshot
	std::vector<particle*> particles;

	// Create a walker from a given particle set 
	walker(std::vector<particle*> particles)
	{
		++ count;
		this->particles = particles;
	}

	// Returns the set of particles in our system
	static std::vector<particle*> create_system_particles()
	{
		// Copy the template system
		std::vector<particle*> ret;
		for (int i=0; i<template_system.size(); ++i)
			ret.push_back(template_system[i]->copy());
		return ret;
	}
};
int walker::count = 0;

// Split a string on whitespace
std::vector<std::string> split_whitespace(std::string to_split)
{
	// Iterator magic
	std::istringstream buffer(to_split);
	std::vector<std::string> split;
	std::copy(std::istream_iterator<std::string>(buffer),
		  std::istream_iterator<std::string>(),
		  std::back_inserter(split));
	return split;
}

// Parse an atom from the input file. Creates a
// nucleus and the appropriate number of electrons
// in the template system.
void parse_atom(std::vector<std::string> split)
{
	double charge = std::stod(split[1]);
	int electrons = std::stoi(split[2]);
	double mass   = std::stod(split[3]);

	// Create the nucleus
	auto n = new nucleus(charge, mass);
	for (int i=4; i<4+dimensions; ++i)
		n->coords[i] = std::stod(split[i]);

	template_system.push_back(n);

	// Add the specified number of electrons
	for (int i=0; i<electrons; ++i)
	{
		auto e = new electron();
		for (int i=0; i<dimensions; ++i)
			e->coords[i] = n->coords[i];
		template_system.push_back(e);
	}
}

// Parse the input file.
void read_input(int np, int pid)
{
	std::ifstream input("input");
	for (std::string line; getline(input, line); )
	{
		auto split = split_whitespace(line);
		std::string tag = split[0];

		// Ignore comments
		if (tag.rfind("!" , 0) == 0) continue;
		if (tag.rfind("#" , 0) == 0) continue;
		if (tag.rfind("//", 0) == 0) continue;

		// Read in the dimensionality
		if (tag == "dimensions")
			dimensions = std::stoi(split[1]);

		// Read in the number of DMC walkers and convert
		// to walkers-per-process
		else if (tag == "walkers") 
			target_population = std::stoi(split[1])/np;

		// Read in the number of DMC iterations
		else if (tag == "iterations")
			dmc_iterations = std::stoi(split[1]);

		// Read in the DMC timestep
		else if (tag == "tau")
			tau = std::stod(split[1]);

		// Read in an atom
		else if (tag == "atom")
			parse_atom(split);

		// Read in an electron
		else if (tag == "electron")
			template_system.push_back(new electron());

		// Adds a noninteracting fermion into the system
		else if (tag == "nif")
			template_system.push_back(new non_interacting_fermion());

		// Add a harmonic well to the system
		else if (tag == "harmonic_well")
			potentials.push_back(new harmonic_well(std::stod(split[1])));
	}

	std::cout <<   "  PID: " << pid 
		  << "\n    Dimensionallity   : " << dimensions
		  << "\n    Iterations        : " << dmc_iterations
	          << "\n    Target population : " << target_population
		  << "\n    Particles         : " << template_system.size()
		  << "\n    Potentials        : " << potentials.size() << "\n";
	input.close();
}

// Returns the number of walkers that should survive after a 
// walker moves from potential pot_before to pot_after. i.e returns
// 0 if the walker should die and 1 if it should live, 2 if another should
// spawn etc...
int walkers_surviving(double pot_before, double pot_after)
{
	double p = exp(-tau*(pot_before + pot_after - 2*trial_energy)/2);
	return std::min(int(floor(p+rand_uniform())), 3);
}

// Contains the results of a single DMC iteration
struct iter_result
{
	int population;
	double average_potential;
	double trial_energy;
	
	// Write the resuls to the "out" file
	void write()
	{
		out_file << population << "," << trial_energy << "\n";
	}

	// Reduce the results of an iteration across processes
	void mpi_reduce(int np, int pid)
	{
		double av_trial_e = 0;
		double av_pot     = 0;
		int pop_sum = 0;
		int ierr    = 0;

		// Sum population over processes
		ierr = MPI_Reduce(&population, &pop_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (ierr != 0) exit(MPI_ERROR);

		// Average average_potential over processes
		ierr = MPI_Reduce(&average_potential, &av_pot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (ierr != 0) exit(MPI_ERROR);

		// Average trial_energy over processes
		ierr = MPI_Reduce(&trial_energy, &av_trial_e, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (ierr != 0) exit(MPI_ERROR);

		// Record results on root process
		if (pid == 0)
		{
			population        = pop_sum;
			average_potential = av_pot     / double(np);
			trial_energy      = av_trial_e / double(np);
		}
	}
};

// The main DMC algorithm
int main(int argc, char** argv)
{
	// Initialize mpi
	if (MPI_Init(&argc, &argv) != 0) exit(MPI_ERROR);

	// For basic timing
	int start_clock = clock();

	// Get the number of processes and my id within them
	int np = 1;
	int pid = 0;
	if (MPI_Comm_size(MPI_COMM_WORLD, &np)  != 0) exit(MPI_ERROR);
	if (MPI_Comm_rank(MPI_COMM_WORLD, &pid) != 0) exit(MPI_ERROR);

	if (pid == 0) 
	{
		std::cout << np << " processes reporting for duty:\n";
		out_file.open("out");
	}

	// Seed random number generator
	srand(pid*clock());

	// Read our input and setup parameters accordingly 
	// do for each process sequentially to avoid access issues
	for (int pid_read = 0; pid_read < np; ++ pid_read)
	{
		if (pid == pid_read) read_input(np, pid);
		MPI_Barrier(MPI_COMM_WORLD);
	}

	// Our DMC walkers
	std::vector<walker*> walkers;

	// Initialize the walkers
	for (int i=0; i<target_population; ++i)
	{
		walker* w = new walker();
		walkers.push_back(w);

		// Carry out an initial large diffusion
		// to avoid unneccasary equilibriation
		// and exact particle overlap at the origin
		w->diffuse(1.0);
	}
	
	// Run our DMC iterations
	for (int iter = 1; iter <= dmc_iterations; ++iter)
	{
		// Output progress on root node
		if (pid == 0) std::cout << "\rIteration: " << iter << "/" << dmc_iterations << "     ";
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
			ir.average_potential += surviving*pot_after;
		}

		// Set our walkers to the newly generated ones
		walkers = new_walkers;

		// Record the resuts of the iteration
		ir.population = walkers.size();
		ir.average_potential /= ir.population;
		ir.trial_energy = trial_energy;

		MPI_Barrier(MPI_COMM_WORLD); // Make sure processes are all on the same iteration
		ir.mpi_reduce(np, pid);      // Accumulate results from all processes
		if (pid == 0) ir.write();    // Output on root process

		// Set trial energy to control population, but allow some fluctuation
		double log_pop_ratio = log(walkers.size()/double(target_population));
		trial_energy = ir.average_potential - log_pop_ratio; 
	}

	// Sample the final wavefunction to file
	if (pid == 0) remove("wavefunction");  // Remove the old wavefunction file
	for (int pid_sample = 0; pid_sample < np; ++ pid_sample)
	{
		if (pid == pid_sample)
		{
			wavefunction_file.open("wavefunction",
				std::ofstream::out | std::ofstream::app);
			for (int n=0; n<walkers.size(); ++n)
				walkers[n]->sample_wavefunction();
			wavefunction_file.close();
		}

		// One process at a time, to avoid access conflicts
		MPI_Barrier(MPI_COMM_WORLD);
	}
	if (pid == 0) std::cout << "\n";

	// Clean up memory
	for (int i=0; i<walkers.size(); ++i) delete walkers[i];
	for (int i=0; i<template_system.size(); ++i) delete template_system[i];
	for (int i=0; i<potentials.size(); ++i) delete potentials[i];
	if  (pid == 0) out_file.close();

	// Output info on objects that werent deconstructed properly
	if (walker::count != 0 || particle::count != 0)
	std::cout << "PID: " << pid << " un-deleted objects:\n"
		  << "  Walkers   : " << walker::count   << "\n"
		  << "  Particles : " << particle::count << "\n";

	// Output success message
	if (pid == 0)
	{
		double time = (clock() - start_clock)/double(CLOCKS_PER_SEC);
		std::cout << "Done, total time: " << time << "s.\n";
	}

	// End mpi safely
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}
