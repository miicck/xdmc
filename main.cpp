#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>

const double PI = 3.14159265358979323846; // Ratio of circumfrence of circle to diameter
const double H_TO_EV = 27.2114;		  // Conversion factor from Hartree to eV
const double EPS_X = 0.001;		  // Distance used for real-space finite differences

double tau = 0.01;		// DMC timestep
int electrons = 1;		// Number of electrons in the system
int target_population = 1000;	// Number of DMC walkers in the simulation
int dmc_iterations = 10000;	// Number of DMC iterations to carry out
double trial_energy = 0;	// Energy used to control population

// Returns a uniform random number in [0,1)
double rand_uniform()
{
	double u = (double)rand();
	return u / RAND_MAX;
}

// Returns a normally distributed random
// number with mean 0 and variance var
double rand_normal(double var)
{
	double u1 = rand_uniform();
	double u2 = rand_uniform();
	return sqrt(-2.0 * var * log(u1))*sin(2.0 * PI * u2);
}

// Represents a single walker in a DMC simulation
class walker
{
public:
	// Default constructor
	walker()
	{
		coords = new double[electrons*3]; // Initialize the coordinates of the electrons
		diffuse(1.0); 			  // Diffuse the electrons to avoid singularities
	}

	// Returns an exact copy of this walker
	walker copy()
	{
		walker ret;

		// Copy my electron coordinates
		for (int i=0; i<electrons*3; ++i)
			ret.coords[i] = coords[i];

		return ret;
	}

	// Diffuse the electrons by moving each of their coordinates
	// by an amount drawn from a normal distribution with variance tau_diff
	void diffuse(double tau_diff)
	{
		for (int i=0; i<electrons*3; ++i)
			coords[i] += rand_normal(tau_diff);
	}

	// Evaluate the potential felt by the electrons in their current locations
	double potential()
	{
		double pot = 0;

		// Evaluate nuclear potential
		for (int i=0; i<electrons; ++i)
		{
			double r = 0;
			for (int d=0; d<3; ++d)
			{
				double xd = coords[i*3+d];
				r += xd*xd;
			}
			r = sqrt(r);
			pot -= 1/r;
		}

		return pot;
	}

	// Evaluate the kinetic energy of this configuration of electrons using
	// the virial theorem, and finite difference derivatives
	// T = 0.5 \sum_i x_i dv/dx_i 
	double kinetic()
	{
		double v0 = potential();
		double kin = 0;

		for (int i=0; i<electrons*3; ++i)
		{
			coords[i] += EPS_X;
			double dv = potential() - v0;
			coords[i] -= EPS_X;
			kin += (coords[i] + EPS_X/2) * dv/EPS_X;
		}
	
		return kin/2;
	}

private:
	double* coords; // The coordinates of the electrons
};

// An object representing the results of a single DMC iteration
struct iter_result
{
	double av_potential; // The average walker potential energy this iteration
	double av_kinetic;   // The average walker kinetic energy this iteration
	int population;      // The population of walkers in this iteration
	int pop_change;      // The population change after this iteration
	int iter_number;     // The iteration number

	// The column titles for the output
	static std::string out_titles()
	{
		std::stringstream ret;
		ret << "Population,Population change,Total energy,Potential energy,Kinetic energy";
		return ret.str();
	}

	// If the column contains data that equilibriates
	static std::string out_equil()
	{
		std::stringstream ret;
		ret << "0,0,1,1,1";
		return ret.str();
		
	}

	// The results of this iteration as a raw string
	std::string out()
	{
		std::stringstream ret;
		ret << population << "," << pop_change << ",";
		ret << av_potential+av_kinetic << "," << av_potential << "," << av_kinetic;
		return ret.str();
	}

	// Print the results of this iteration nicely formatted
	void print()
	{
		std::cout << "\nIteration " << iter_number + 1 << "/" << dmc_iterations << "\n";
		std::cout << "Population   : " << population << "\n";
		std::cout << "Energy       : " << av_potential + av_kinetic << "\n";
		std::cout << "    Potential: " << av_potential << "\n";
		std::cout << "    Kinetic  : " << av_kinetic << "\n";
	}
};

std::vector<walker> walkers(target_population); // The DMC walkers
std::vector<iter_result> iter_results;		// The results of DMC iterations accumulated so far

// Returns the number of walkers that survive after a walker
// moves from a configuration with potential pot_before
// to a configuration with potential pot_after. A result
// of 0 indicates this walker should die, a result of 1
// that it should live and a result of 2 that it should
// live and another walker be spawned in the same place.
int walkers_surviving(double pot_before, double pot_after)
{
	double av_pot = (pot_before + pot_after)/2;
	double p = exp(-tau * (av_pot - trial_energy));

	if (p < 1)
	{
		if (rand_uniform() < p) return 1;
		return 0;
	}
	if (rand_uniform() < p-1) return 2;
	return 1;
}

// Carry out a single DMC iteration
void iterDMC()
{
	// Initialize the results of this iteration
	iter_result ir;
	ir.av_potential = 0;
	ir.av_kinetic   = 0;
	ir.population = walkers.size();

	// Will contain new walkers that were spawned this iteration
	std::vector<walker> walkers_spawned;

	for (int n=walkers.size()-1; n>=0; --n)
	{
		// Evaluate the potential and kinetic energy
		// before the diffusion step
		double pot_before = walkers[n].potential();
		double kin_before = walkers[n].kinetic();

		// Perform a diffusion step
		walkers[n].diffuse(tau);

		// Evaluate the potential and kinetic energy
		// after the diffusion step
		double pot_after  = walkers[n].potential();
		double kin_after = walkers[n].kinetic();

		// Sum averages
		ir.av_kinetic += (kin_before+kin_after)/2;
		ir.av_potential += (pot_before+pot_after)/2;

		// Work out how many walkers survive the move
		int surviving = walkers_surviving(pot_before, pot_after);
		if (surviving == 0) walkers.erase(walkers.begin()+n); // Kill this walker
		else if (surviving > 1) walkers_spawned.push_back(walkers[n].copy()); // Spawn a new walker
	}	

	// Append the walkers spawned to the collection of all walkers
	walkers.insert(walkers.end(), walkers_spawned.begin(), walkers_spawned.end());

	// Evaluate averages
	ir.av_potential /= ir.population; 
	ir.av_kinetic   /= ir.population;
	ir.iter_number   = iter_results.size();
	ir.pop_change    = walkers.size() - ir.population;
	iter_results.push_back(ir);

	// Set the trial energy to the average potential
	// so that the population changes as little as possible
	trial_energy = ir.av_potential;

	// Print the results of the iteration
	ir.print();
}

// Output the results of the DMC calculation to disk 
void outputDMC()
{
	// Accumulate the iteraion results
	std::stringstream to_write;
	for (int n=0,ntest=iter_results.size(); n<ntest; ++n)
		to_write << iter_results[n].out() << "\n";
 
	// Write them to disk
	std::ofstream output("out");
	std::string titles = iter_result::out_titles();
	std::string equils = iter_result::out_equil();
	titles.append("\n");
	equils.append("\n");
	output.write(titles.c_str(),titles.size());
	output.write(equils.c_str(),equils.size());
	output.write(to_write.str().c_str(),to_write.str().size());
}

// Program entrypoint
int main(int argc, char** argv)
{
	srand(clock());
	for (int n=0; n<dmc_iterations; ++n)
		iterDMC();
	outputDMC();
}
