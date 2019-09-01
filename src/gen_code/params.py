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
},{
    "in_name"     : "walkers",
    "type"        : "unsigned",
    "cpp_name"    : "target_population",
    "default"     : "1000",
    "description" : ("The target population of DMC walkers. The actual number of "
                     "dmc walkers will fluctuate during runtime, but will be bias "
                     "towards this value."),
},{
    "in_name"     : "iterations",
    "type"        : "unsigned",
    "cpp_name"    : "dmc_iterations",
    "default"     : "10000",
    "description" : ("The number of DMC iterations, each corresponding to a "
                     "step of tau in imaginary time."),
},{
    "in_name"     : "iter",
    "type"        : "unsigned",
    "cpp_name"    : "dmc_iteration",
    "default"     : "0",
    "description" : "The current DMC iteration.",
},{
    "in_name"     : "tau",
    "type"        : "double",
    "cpp_name"    : "tau",
    "default"     : "0.01",
    "description" : "The DMC timestep in atomic units.",
},{
    "in_name"     : "tau_psi",
    "type"        : "double",
    "cpp_name"    : "tau_psi",
    "default"     : "0.1",
    "description" : ("The effective timestep used when evaluating the diffused "
                     "wavefunction. Larger values encourage the formation of nodal "
                     "pockets.")
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
    "description" : ("The DMC trial energy in atomic units (Hartree). This "
                     "value is used to control the DMC population and will "
                     "fluctuate during runtime. After equilibriation, it will "
                     "fluctuate around the ground state energy of the system."),
},{
    "in_name"     : "pre_diffusion",
    "type"        : "double",
    "cpp_name"    : "pre_diffusion",
    "default"     : "1.0",
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
    "default"     : "true",
    "description" : "True if nodal surface files are to be written.",
},{
    "in_name"     : "exact_1d_nodes",
    "type"        : "bool",
    "cpp_name"    : "exact_1d_nodes",
    "default"     : "false",
    "description" : ("True if the exact nodal surface is used when we're in 1 "
                     "spatial dimension.")
},{
    "in_name"     : "exchange_prob",
    "type"        : "double",
    "cpp_name"    : "exchange_prob",
    "default"     : "0.5",
    "description" : ("The probability of a walker making an exchange move in any "
                     "given timestep. The actual exchange move made will be chosen "
                     "at random. 1 - this is the probability of simply diffusing, "
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
}]

