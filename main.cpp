#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>

// Forward declerations
class walker;
class permutor;
struct iter_result;

const double PI = 3.14159265358979323846; // Ratio of circumfrence of circle to diameter
const double H_TO_EV = 27.2114;		  // Conversion factor from Hartree to eV
const double EPS_X = 0.001;		  // Distance used for real-space finite differences

std::vector<walker> walkers; 		// The DMC walkers
std::vector<iter_result> iter_results;	// The results of DMC iterations accumulated so far

double tau = 0.01;		// DMC timestep
int target_population = 1000;	// Number of DMC walkers in the simulation
int dmc_iterations = 1000;	// Number of DMC iterations to carry out
double trial_energy = 0;	// Energy used to control population
int electrons = 1;		// Number of electrons in the system
int nuclei    = 1;		// The number of nuclei in the system
double* nuclear_coords;		// The coordinates of nucleii
double* nuclear_charge;		// The nuclear charges
permutor* permutations;		// Will store permutations of 1,...,electrons
bool use_permutations = true;   // Use electron coordinate permutations?

// Returns a uniform random number in [0,1)
double rand_uniform()
{
	double u = (double)rand();
	return u / RAND_MAX;
}

// Returns a normally distributed random
// number with mean 0 and variance var
// uses a box-muller transform
double rand_normal(double var)
{
	double u1 = rand_uniform();
	double u2 = rand_uniform();
	return sqrt(-2.0 * var * log(u1))*sin(2.0 * PI * u2);
}

// Computes the factorial of n
int factorial(int n)
{
	if (n==0) return 1;
	return n*factorial(n-1);
}

// Computes all of the permutations of 1,...,n and stores the
// result in permutations. Additionaly the sign of each permutation
// is stored. Will allocate memory for the permutations and 
// signs variables. Uses heaps algorithm.
class permutor
{
public: 
	int** perms;
	int* signs;
	int fact;
	int n;

	void print()
	{
		for (int i=0; i<fact; ++i)
		{
			for (int j=0; j<n; ++j)
				std::cout << perms[i][j];
			std::cout << "\n";
		}
	}

	permutor(int number)
	{
		n = number;
		fact = factorial(n);
		int current_row = 1;

		int* a = new int[n];
		int* c = new int[n];

		signs = new int[fact];
		perms = new int*[fact];
		for (int i=0; i<fact; ++i)
			perms[i] = new int[n];


		for (int i=0; i<n; ++i)
		{
			a[i] = i;
			c[i] = 0;
			perms[0][i] = i;
			signs[0]    = 1;
		}

		int i=0;
		int sgn=1;
		while(i < n)
		{
			if (c[i] < i)
			{
				if (i%2 == 0)
				{
					int temp = a[i];
					a[i] = a[0];
					a[0] = temp;
					sgn  = -sgn;
				}
				else
				{
					int temp = a[i];
					a[i]     = a[c[i]];
					a[c[i]]  = temp;
					sgn      = -sgn;
				}

				for (int m=0; m<n; ++m)
					perms[current_row][m] = a[m];
				signs[current_row] = sgn;
				++ current_row;
				++ c[i];
				i = 0;
			}
			else
			{
				c[i] = 0;
				++ i;
			}
		}
	}
};

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

		int pmax = 1;
		if (use_permutations) pmax = permutations->fact;
		for (int p=0; p<pmax; ++p)
		{
			double pot_conf = 0;

			// Electron-nuclear interactions
			for (int i=0; i<electrons; ++i)
				for (int j=0; j<nuclei; ++j)
				{
					double r = 0;
					for (int d=0; d<3; ++d)
					{
						int ip = i;
						if (use_permutations) ip = permutations->perms[p][i];
						double xd = coords[ip*3+d] - nuclear_coords[j*3+d];
						r += xd*xd;
					}
					r = sqrt(r);
					pot_conf -= nuclear_charge[j]/r;
				}

			// Electron-electron interactions
			for (int i=0; i<electrons; ++i)
				for (int j=0; j<i; ++j)
				{
					double r = 0;
					for (int d=0; d<3; ++d)
					{
						int ip = i; 
						int jp = j;
						if (use_permutations)
						{
							ip = permutations->perms[p][i];
							jp = permutations->perms[p][j];
						}
						double xd = coords[ip*3+d] - coords[jp*3+d];
						r += xd*xd;
					}
					r = sqrt(r);
					pot_conf += 1/r;
				}

			if (use_permutations) pot += pot_conf * permutations->signs[p];
			else pot += pot_conf;
		}

		return pot/pmax;
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
	double log_pop_ratio = log(double(target_population)/double(walkers.size()));
	trial_energy = ir.av_potential + log_pop_ratio;
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

// Read the input file, returns true if successful
bool read_input(char* filename)
{
	std::ifstream file(filename);
	if (file.fail())
	{
		std::cout << "Could not open input file!\n";
		return false;
	}

	// Parse the input file line by line
	// each line will either be of the form
	// "keyword keyword_value"
	// or
	// "Z x y z"
	// where Z is the charge of an atom at
	// x,y,z that we wish to include in the
	// simulation
	std::cout << "Input file:" << "\n";
	std::string line;
	std::vector<double> input_doubles; // Will contain the nuclear information
	while(std::getline(file, line))
	{
		std::cout << "    " << line << "\n";
		std::stringstream ss(line);
		std::string word;
		while(std::getline(ss,word,' '))
		{
			if (word == "walkers")
			{
				// Read in the population
				std::getline(ss,word,' ');
				target_population = std::stoi(word);
			}
			else if (word == "iterations")
			{
				// Read in the number of DMC iterations
				std::getline(ss,word,' ');
				dmc_iterations = std::stoi(word);
			}
			else if (word == "timestep")
			{
				// Read in the DMC timestep
				std::getline(ss,word,' ');
				tau = std::stod(word);
			}
			else if (word == "use_permutations")
			{
				// Read in if we use permutations
				std::getline(ss,word,' ');
				std::istringstream(word.c_str()) >> use_permutations;
			}
			else input_doubles.push_back(std::stod(word));	
		}
	}
	std::cout << "End input file\n\n";
	
	// Create the nuclei
	nuclei = input_doubles.size()/4;
	nuclear_coords = new double[nuclei*3];
	nuclear_charge = new double[nuclei];
	double total_charge = 0;
	for (int i=0; i<nuclei; ++i)
	{
		nuclear_charge[i] = input_doubles[i*4];
		total_charge += nuclear_charge[i];
		for (int d=0; d<3; ++d)
			nuclear_coords[i+d] = input_doubles[i*4+d+1];
	}

	// Make the system neutral
	electrons = (int)total_charge;
	permutations = new permutor(electrons);

	// Create the walkers representing the electrons
	for (int i=0; i<target_population; ++i)
	{
		walker w;
		walkers.push_back(w);
	}
	
	// Seed the random number generator
	int rand_seed = clock();
	srand(rand_seed);
	
	// Output the input paremeters
	std::cout << "Input parameters:\n";
	std::cout << "    Nuclei             : " << nuclei << "\n";
	std::cout << "    Total charge       : " << total_charge << "\n"; 
	std::cout << "    Electrons          : " << electrons << "\n";
	std::cout << "    Walkers            : " << target_population << "\n";
	std::cout << "    Iterations         : " << dmc_iterations << "\n";
	std::cout << "    DMC timestep       : " << tau << "\n";
	std::cout << "    Random seed        : " << rand_seed << "\n";
	std::cout << "    Using permutations : " << use_permutations << "\n"; 
	std::cout << "End input parmeters\n\n";

	return true;
}

// Program entrypoint
int main(int argc, char** argv)
{
	if (argc < 1)
	{
		std::cout << "No input file given!\n";
		return -1;
	}

	if (!read_input(argv[1])) return -2;

	for (int n=0; n<dmc_iterations; ++n)
	{
		iterDMC();
		std::cout << "\rRunning dmc, progress: " << n+1 <<"/" << dmc_iterations;
		std::cout << " (" << 100*double(n+1)/double(dmc_iterations) << "%)";
		iter_result last = iter_results[iter_results.size()-1];
		std::cout << ". Last energy: " << last.av_potential + last.av_kinetic;
		std::cout << " Last population: " << last.population;
		std::cout << "               ";
	}
	outputDMC();
}
