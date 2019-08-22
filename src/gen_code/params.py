# This file acts as a central source for all of the parameter definitions. It is
# used to generate the parameters in the cpp code, generate documentation and to
# generate automatic parameter parsing/output code.
#
# The syntax for a parameter is
#
#   params["name_read_from_input_file"] = [
#     "cpp_datatype",               # String - directly inserted into c++ file(s)
#     "cpp_variable_name",          # String - directly inserted into c++ file(s)
#     "default_value",              # String - directly inserted into c++ file(s)
#     "description",                # String - directly inserted into c++ file(s)
#      is_read_from_input_file      # Bool   - used to inform code generation
#   ]
#

params = {}

params["pid"] = [
    "int",
    "pid",
    "0",
    "The MPI process id of this process. Will be in [0,np).",
    False
]

params["np"] = [
    "int",
    "np",
    "1",
    "The number of MPI processes.",
    False
]

params["dimensions"] = [
    "int",
    "dimensions",
    "3",
    "The spatial dimensions of the system.",
    True
]

params["walkers"] = [
    "int",
    "target_population",
    "1000",
    ("The target population of DMC walkers. The actual number of "
     "dmc walkers will fluctuate during runtime, but will be bias "
     "towards this value."),
     True
]

params["max_pop_ratio"] = [
    "double",
    "max_pop_ratio",
    "4.0",
    ("The maximum allowed population, expressed as a multiple of "
     "the target population."),
     True
]

params["min_pop_ratio"] = [
    "double",
    "min_pop_ratio",
    "0.0",
    ("The minimum allowed population, expressed as a fraction of "
    "the target population."),
    True
]

params["iterations"] = [
    "int",
    "dmc_iterations",
    "10000",
    ("The number of DMC iterations, each corresponding to a "
     "step of tau in imaginary time."),
     True
]

params["iter"] = [
    "int",
    "dmc_iteration",
    "0",
    ("The current DMC iteration."),
    False
]

params["tau"] = [
    "double",
    "tau",
    "0.01",
    "The DMC timestep in atomic units.",
    True
]

params["tau_c_ratio"] = [
    "double",
    "tau_c_ratio",
    "1.0",
    ("The ratio of tau_c:tau, where tau_c is the effective "
     "cancellation timestep and tau is the DMC timestep."),
     True
]

params["trial_energy"] = [
    "double",
    "trial_energy",
    "0.0",
    ("The DMC trial energy in atomic units (Hartree). This "
     "value is used to control the DMC population and will "
     "fluctuate during runtime. After equilibriation, it will "
     "fluctuate around the ground state energy of the system."),
     True
]

params["best_energy"] = [
    "double",
    "best_energy",
    "0.0",
    ("The best estimate of the DMC energy in atomic units (Hartree). "
     "This will be updated as the simulation progresses."),
     True
]

params["pre_diffusion"] = [
    "double",
    "pre_diffusion",
    "1.0",
    ("The amount of imaginary time that the walkers will diffuse for "
     "before the first full DMC iteration. Effectively, this is how "
     "spread out the initial wavefunction is."),
     True
]  

params["write_wavefunction"] = [
   "bool",
   "write_wavefunction",
   "true",
   "True if wavefunction files are to be written.",
   True
]

params["exchange_moves"] = [
    "bool",
    "exchange_moves",
    "true",
    "True if exchange moves are to be made.",
    True
]

params["exchange_prob"] = [
    "double",
    "exchange_prob",
    "0.5",
    ("The probability of a walker making an exchange move in any "
     "given timestep. The actual exchange move made will be chosen at random. "
     "1 - this is the probability of simply diffusing, making no exchange moves."),
     True
]

params["cancel_scheme"] = [
    "std::string",
    "cancel_scheme",
    '"diffusive"',
    "The cancellation scheme used.",
    True
]

params["cancel_function"] = [
    "std::string",
    "cancel_function",
    '"integrated_weight"',
    "The cancellation function used if we are applying pairwise cancellations.",
    True
]

params["cancelled_weight"] = [
    "double",
    "cancelled_weight",
    "0.0",
    ("The amount of weight cancelled during the last cancellation step."),
    False
]

params["correct_seperations"] = [
    "bool",
    "correct_seperations",
    "false",
    "True if seperation corrections are applied.",
    True
]
