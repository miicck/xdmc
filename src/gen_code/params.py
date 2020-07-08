# 
#     XDMC
#     Copyright (C) Michael Hutcheon (email mjh261@cam.ac.uk)
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     For a copy of the GNU General Public License see <https://www.gnu.org/licenses/>.
# 
params = [{
    "type"        : "std::string",
    "cpp_name"    : "run_mode",
    "default"     : '"calculation"',
    "description" : "The type of execution running.",
},{
    "in_name"     : "dimensions",
    "type"        : "unsigned",
    "cpp_name"    : "dimensions",
    "default"     : "3",
    "allowed"     : "positive",
    "description" : "Spatial dimensions of the system (i.e 1D, 2D, 3D etc..).",
},{
    "in_name"     : "walkers",
    "type"        : "unsigned",
    "cpp_name"    : "target_population",
    "default"     : "1000",
    "allowed"     : "positive",
    "description" : ("The target population of DMC walkers. The actual number of "
                     "dmc walkers will fluctuate during runtime, but will be bias "
                     "towards this value."),
},{
    "in_name"     : "walkers_end",
    "type"        : "unsigned",
    "cpp_name"    : "end_target_population",
    "default"     : "0",
    "allowed"     : "positive",
    "description" : ("The target number of walkers at the end of the simulation. "
                     "The target will decrease linearly from the initial value "
                     "specified by the walkers keyword to this value. A value of "
                     "0 will lead to a constant population. ")
},{
    "in_name"     : "iterations",
    "type"        : "unsigned",
    "cpp_name"    : "dmc_iterations",
    "default"     : "10000",
    "allowed"     : "positive",
    "description" : ("The number of DMC iterations, each corresponding to a "
                     "step of tau in imaginary time."),
},{
    "type"        : "unsigned",
    "cpp_name"    : "dmc_iteration",
    "default"     : "0",
    "description" : "The current DMC iteration.",
},{
    "in_name"     : "tau",
    "type"        : "double",
    "cpp_name"    : "tau",
    "default"     : "0.01",
    "allowed"     : "positive",
    "description" : "The DMC timestep in atomic units.",
},{
    "in_name"     : "tau_nodes",
    "type"        : "double",
    "cpp_name"    : "tau_nodes",
    "default"     : "params::tau",
    "allowed"     : "positive",
    "description" : ("The effective timestep used enforce fermionic antisymmetry. "
                     "Essentially corresponds to the range of influence of a single "
                     "walker. Setting tau_nodes >> tau will combat bosonic collapse  "
                     "(but may introduce bias).")
},{
    "in_name"     : "tau_nodes_estimator",
    "type"        : "std::string",
    "cpp_name"    : "tau_nodes_estimator",
    "default"     : '"none"',
    "allowed"     : "strings none min_sep min_sep_mpi",
    "description" : "Method used to estimate tau_nodes during runtime."
},{
    "in_name"     : "self_gf_strength",
    "type"        : "double",
    "cpp_name"    : "self_gf_strength",
    "default"     : "1.0",
    "allowed"     : "between 0.0 1.0",
    "description" : ("How much a walker contributes to it's own diffused "
                     "wavefunction. 1.0 <=> the same as other walkers. 0.0 "
                     " <=> not at all.")
},{
    "in_name"     : "energy_estimator",
    "type"        : "std::string",
    "cpp_name"    : "energy_estimator",
    "default"     : '"growth"',
    "allowed"     : "strings growth potential",
    "description" : "The method used to estimate the energy."
},{
    "in_name"     : "growth_mixing_factor",
    "type"        : "double",
    "cpp_name"    : "growth_mixing_factor",
    "default"     : "0.0",
    "allowed"     : "between 0.0 1.0",
    "description" : ("The amount that the old trial energy is mixed back "
                     "into the new trial energy each iteration.")
},{
    "in_name"     : "full_exchange",
    "type"        : "bool",
    "cpp_name"    : "full_exchange",
    "default"     : "false",
    "description" : ("If true, sample exchange moves from permutations, "
                     "rather than permutation operators.")
},{
    "in_name"     : "correct_average_weight",
    "type"        : "bool",
    "cpp_name"    : "correct_average_weight",
    "default"     : "false",
    "description" : ("If true, additional exchange moves will be made in an "
                     "attempt to keep the average walker weight = 0.")
},{
    "in_name"     : "coulomb_softening",
    "type"        : "double",
    "cpp_name"    : "coulomb_softening",
    "default"     : "0",
    "allowed"     : "positive",
    "description" : ("Softening parameter for the coulomb potential. 0 corresponds "
                     "to a bare coulomb potential, larger values correspond to softer "
                     "potentials of the form 1/(r + s).")
},{
    "in_name"     : "diffusion_scheme",
    "type"        : "std::string",
    "cpp_name"    : "diffusion_scheme",
    "default"     : '"stochastic_nodes_mpi"',
    "allowed"     : ("strings exact_1d bosonic exchange_diffuse "
                     "max_seperation max_seperation_mpi "
                     "stochastic_nodes stochastic_nodes_mpi "
                     "stochastic_nodes_permutations"),
    "description" : "The diffusion scheme used."
},{
    "in_name"     : "max_weight",
    "type"        : "double",
    "cpp_name"    : "max_weight",
    "default"     : "10",
    "allowed"     : "positive",
    "description" : ("The maximum allowed weight of any individual walker. "
                     "If this is exceeded then the iteration is reverted. Smaller "
                     "values combat population explosion, but introduce bias. Very "
                     "small values may lead to excessive numbers of reverted "
                     "iterations. If working with a system with coulomb interactions "
                     "use max_weight in conjunction with coulomb_softening for best "
                     "results.")
},{
    "type"        : "int",
    "cpp_name"    : "np",
    "default"     : "1",
    "description" : "The number of MPI processes.",
},{
    "type"        : "int",
    "cpp_name"    : "pid",
    "default"     : "0",
    "description" : "The MPI process id of this process. Will be in [0,np).",
},{
    "in_name"     : "trial_energy",
    "type"        : "double",
    "cpp_name"    : "trial_energy",
    "default"     : "0.0",
    "description" : ("The initial value of the DMC trial energy in atomic units (Hartree). "
                     "This value is used to control the DMC population and will "
                     "fluctuate during runtime. After equilibriation, it will "
                     "fluctuate around the ground state energy of the system."),
},{
    "type"        : "double",
    "cpp_name"    : "cancelled_weight",
    "default"     : "0.0",
    "description" : ("The total amount of weight cancelled in the last iteration.")
},{
    "in_name"     : "pre_diffusion",
    "type"        : "double",
    "cpp_name"    : "pre_diffusion",
    "default"     : "1.0",
    "alowed"      : "positive",
    "description" : ("The amount of imaginary time that the walkers will diffuse for "
                     "before the first full DMC iteration. Effectively, this is how "
                     "spread out the initial wavefunction is."),
},{
    "in_name"     : "write_wavefunction",
    "type"        : "bool",
    "cpp_name"    : "write_wavefunction",
    "default"     : "true",
    "description" : "True if wavefunction files are to be written.",
},{
    "in_name"     : "write_nodal_surface",
    "type"        : "bool",
    "cpp_name"    : "write_nodal_surface",
    "default"     : "false",
    "description" : "True if nodal surface files are to be written.",
},{
    "in_name"     : "exchange_prob",
    "type"        : "double",
    "cpp_name"    : "exchange_prob",
    "default"     : "0.5",
    "allowed"     : "between 0.0 1.0",
    "description" : ("The probability of a walker making an exchange move in any "
                     "given timestep. The actual exchange move made will be chosen "
                     "at random. (1 - this) is the probability of simply diffusing, "
                     "making no exchange moves."),
},{
    "type"        : "int",
    "cpp_name"    : "start_clock",
    "default"     : "0",
    "description" : "The result of clock() called at program start."
},{
    "type"        : "int",
    "cpp_name"    : "dmc_start_clock",
    "default"     : "0",
    "description" : "The result of clock() called just before first DMC iteration."
},{
    "type"        : "int",
    "cpp_name"    : "argc",
    "default"     : "0",
    "description" : "The number of arguments passed to the program."
},{
    "type"        : "char**",
    "cpp_name"    : "argv",
    "default"     : "nullptr",
    "description" : "The arguments passed to the program."
}]

