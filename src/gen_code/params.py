# This file acts as a central source for all of the parameter definitions. It is
# used to generate the parameters in the cpp code, generate documentation and to
# generate automatic parameter parsing/output code.
# The order of parameter definitions here will be the order they are read in and
# output.

params = [{
    "in_name"     : "dimensions",
    "type"        : "unsigned",
    "cpp_name"    : "dimensions",
    "default"     : "3",
    "description" : "The spatial dimensions of the system.",
    "read_in"     : True
},{
    "in_name"     : "walkers",
    "type"        : "unsigned",
    "cpp_name"    : "target_population",
    "default"     : "1000",
    "description" : ("The target population of DMC walkers. The actual number of "
                     "dmc walkers will fluctuate during runtime, but will be bias "
                     "towards this value."),
    "read_in"     : True
},{
    "in_name"     : "max_pop_ratio",
    "type"        : "double",
    "cpp_name"    : "max_pop_ratio",
    "default"     : "4.0",
    "description" : ("The maximum allowed population, expressed as a multiple of "
                     "the target population."),
    "read_in"     : True
},{
    "in_name"     : "min_pop_ratio",
    "type"        : "double",
    "cpp_name"    : "min_pop_ratio",
    "default"     : "0.0",
    "description" : ("The minimum allowed population, expressed as a fraction of "
                     "the target population."),
    "read_in"     : True
},{
    "in_name"     : "iterations",
    "type"        : "unsigned",
    "cpp_name"    : "dmc_iterations",
    "default"     : "10000",
    "description" : ("The number of DMC iterations, each corresponding to a "
                     "step of tau in imaginary time."),
    "read_in"     : True
},{
    "in_name"     : "iter",
    "type"        : "unsigned",
    "cpp_name"    : "dmc_iteration",
    "default"     : "0",
    "description" : "The current DMC iteration.",
    "read_in"     : False
},{
    "in_name"     : "tau",
    "type"        : "double",
    "cpp_name"    : "tau",
    "default"     : "0.01",
    "description" : "The DMC timestep in atomic units.",
    "read_in"     : True
},{
    "in_name"     : "tau_c_ratio",
    "type"        : "double",
    "cpp_name"    : "tau_c_ratio",
    "default"     : "1.0",
    "description" : ("The ratio of tau_c:tau, where tau_c is the effective "
                     "cancellation timestep and tau is the DMC timestep."),
    "read_in"     : True
},{
    "in_name"     : "np",
    "type"        : "int",
    "cpp_name"    : "np",
    "default"     : "1",
    "description" : "The number of MPI processes.",
    "read_in"     : False
},{
    "in_name"     : "pid",
    "type"        : "int",
    "cpp_name"    : "pid",
    "default"     : "0",
    "description" : "The MPI process id of this process. Will be in [0,np).",
    "read_in"     : False
},{
    "in_name"     : "trial_energy",
    "type"        : "double",
    "cpp_name"    : "trial_energy",
    "default"     : "0.0",
    "description" : ("The DMC trial energy in atomic units (Hartree). This "
                     "value is used to control the DMC population and will "
                     "fluctuate during runtime. After equilibriation, it will "
                     "fluctuate around the ground state energy of the system."),
     "read_in"    : True
},{
    "in_name"     : "pre_diffusion",
    "type"        : "double",
    "cpp_name"    : "pre_diffusion",
    "default"     : "1.0",
    "description" : ("The amount of imaginary time that the walkers will diffuse for "
                     "before the first full DMC iteration. Effectively, this is how "
                     "spread out the initial wavefunction is."),
    "read_in"     : True
},{
    "in_name"     : "write_wavefunction",
    "type"        : "bool",
    "cpp_name"    : "write_wavefunction",
    "default"     : "true",
    "description" : "True if wavefunction files are to be written.",
    "read_in"     : True
},{
    "in_name"     : "exchange_moves",
    "type"        : "bool",
    "cpp_name"    : "exchange_moves",
    "default"     : "true",
    "description" : "True if exchange moves are to be made.",
    "read_in"     : True
},{
    "in_name"     : "exchange_prob",
    "type"        : "double",
    "cpp_name"    : "exchange_prob",
    "default"     : "0.5",
    "description" : ("The probability of a walker making an exchange move in any "
                     "given timestep. The actual exchange move made will be chosen "
                     "at random. 1 - this is the probability of simply diffusing, "
                     "making no exchange moves."),
    "read_in"     : True
},{
    "in_name"     : "cancel_scheme",
    "type"        : "std::string",
    "cpp_name"    : "cancel_scheme",
    "default"     : '"diffusive"',
    "description" : "The cancellation scheme used.",
    "read_in"     : True
},{
    "in_name"     : "cancel_function",
    "type"        : "std::string",
    "cpp_name"    : "cancel_function",
    "default"     : '"integrated_weight"',
    "description" : "The cancellation function used if we are applying pairwise cancellations.",
    "read_in"     : True
},{
    "in_name"     : "cancelled_weight",
    "type"        : "double",
    "cpp_name"    : "cancelled_weight",
    "default"     : "0.0",
    "description" : "The amount of weight cancelled during the last cancellation step.",
    "read_in"     : False
},{
    "in_name"     : "correct_seperations",
    "type"        : "bool",
    "cpp_name"    : "correct_seperations",
    "default"     : "false",
    "description" : "True if seperation corrections are applied.",
    "read_in"     : True
},{
    "in_name"     : "self_sign_strength",
    "type"        : "double",
    "cpp_name"    : "self_sign_strength",
    "default"     : "0.0",
    "description" : ("The strength of a walkers own greens function in "
                     "deciding what sign it should be. A value of 0 strongly "
                     "encourages the formation of nodal pockets, but represents "
                     "an approximation to the true greens function. A value of " 
                     "1 reproduces the correct greens function but is likely to "
                     "lead to nodal pockets breaking apart."),
    "read_in"     : True
}]
