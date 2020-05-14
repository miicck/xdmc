# 
#     XDMC
#     Copyright (C) 2019 Michael Hutcheon (email mjh261@cam.ac.uk)
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
import sys
from parser import parse_evolution
import numpy as np

def exp_fit():
    from scipy.optimize import curve_fit
    for name, energies in zip(*parse_evolution()):
        if not "Trial energy" in name: continue

        # Fit energy vs iteration to exponential decay
        # using initial guess informed by data
        to_fit = lambda n, a, b, c: a + b*np.exp(-n/c)
        p0 = [energies[-1], energies[0] - energies[-1], len(energies)/4]
        bounds = [[min(energies), -np.inf, 1],
                  [max(energies),  np.inf, 10*len(energies)]]
        
        x = np.array([float(i) for i in range(0, len(energies))])
        par, cov = curve_fit(to_fit, x, energies, p0, bounds=bounds)
        residuals = [e - to_fit(n, *par) for n, e in enumerate(energies)]
        print(par[0] ,"+/-", np.std(residuals))

        if "plot" in sys.argv:
            import matplotlib.pyplot as plt
            
            # Plot the data, fit and residuals
            plt.subplot(211)
            plt.plot(energies, label="DMC")
            xs = np.linspace(0, len(energies), 1000)
            ys = [to_fit(x, *par) for x in xs]
            plt.plot(xs, ys, label="Fit")
            plt.ylabel("Trial energy (Hartree)")
            plt.legend()

            plt.subplot(212)
            plt.plot(residuals)
            plt.axhline(0, color="black")
            plt.ylabel("Residuals (Hartree)")

            plt.show()

def multi_exp_fit(exponents=3):
    from scipy.optimize import curve_fit
    for name, energies in zip(*parse_evolution()):
        if not "Trial energy" in name: continue
        to_fit = lambda n, a, b, c: a + b*np.exp(-n/c)

        residuals  = np.array(energies)
        parameters = []
        for i in range(0, exponents):
            p0 = [residuals[-1], residuals[0] - residuals[-1], len(residuals)/4]
            bounds = [[min(residuals), -np.inf, 1],
                      [max(residuals), np.inf,  10*len(residuals)]]
            par, cov = curve_fit(to_fit, range(0, len(residuals)),
                                 residuals, p0, bounds=bounds)
            fit = [to_fit(n, *par) for n in range(0, len(residuals))]
            residuals -= fit
            parameters.append(par)

        model = lambda n : sum(to_fit(n, *p) for p in parameters)
        print(model(10e10))

        if "plot" in sys.argv:
            import matplotlib.pyplot as plt

            fit = [model(n) for n in range(0, len(energies))]
            ext = [model(n) for n in range(0, len(energies)*10)]
            plt.subplot(211)
            plt.plot(energies)
            plt.plot(fit)
            plt.plot(ext)
            plt.subplot(212)
            plt.plot(np.array(energies)-np.array(fit))
            plt.axhline(0, color="black")
            plt.show()

def average():
    for name, energies in zip(*parse_evolution()):
        if not "Trial energy" in name: continue
        start = int(sys.argv[2])

        energies = energies[start:]
        print(np.mean(energies), " +/- ", np.std(energies)/np.sqrt(len(energies)))
        if "plot" in sys.argv:
            import matplotlib.pyplot as plt

            plt.plot(energies)
            plt.axhline(np.mean(energies), color="black")
            plt.show()

def reblocked():
    es = None
    for name, energies in zip(*parse_evolution()):
        if not "Trial energy" in name: continue
        es = energies
    if es is None: raise Exception("Could not find the trial energies!")
    start = int(sys.argv[2])

    # Get the energies after the provided 
    # start iteration, and their mean
    es = es[start:]
    energy = np.mean(es)

    # Init arrays
    block_errors = []
    block_error_errors = []

    # Loop over blocking transformation number
    blocking_transform = 0
    while True:

        # Work out the block length, if it is greater than
        # the total length, break
        block_length = 2 ** blocking_transform 
        blocking_transform += 1
        if block_length > len(es):
            break

        # Work out the energies of each block
        block_energies = []
        offset = 0
        while offset < len(es):
            block_energies.append(np.mean(es[offset:offset+block_length]))
            offset += block_length

        # Work out the reblocked error + the error on that
        block_err = np.std(block_energies) / np.sqrt(len(block_energies))
        block_errors.append(block_err)
        block_err_err = block_err / np.sqrt(2 * (len(block_energies) - 1))
        block_error_errors.append(block_err_err)

    # Plot vs blocking number if requested
    if "plot" in sys.argv:
        import matplotlib.pyplot as plt

        plt.errorbar(range(len(block_errors)), block_errors, yerr=block_error_errors)
        plt.xlabel("Blocking transformation number")
        plt.ylabel("Reblocked errror (Ha)")
        plt.show()

    print("{0} +/- {1}".format(energy, max(block_errors)))

def pi_weighted():
    for name, energies in zip(*parse_evolution()):
        if not "Trial energy" in name: continue
        start = int(sys.argv[2])
        corr  = int(sys.argv[3])
        energies = energies[start:]
        
        f = open("progress")
        lines = f.read().split("\n")
        f.close()

        for l in lines:
                if "DMC timestep" in l:
                        tau = float(l.split(":")[-1])
                        break

        mean_e = np.mean(energies)
        pi = lambda m : np.product([np.exp(-tau*(e-mean_e)) for e in energies[m-corr:m]])
        estimate  = sum([pi(m)*energies[m] for m in range(0, len(energies))])
        estimate /= sum([pi(m) for m in range(0, len(energies))])
        print(estimate)

def pop_correction():
    for name, data in zip(*parse_evolution()):
        if "Population"   in name: pop = np.array(data)
        if "Trial energy" in name: energy = np.array(data)

    half   = int(len(energy)/2)
    energy = energy[half:]
    pop    = pop[half:]
    pop_f  = pop - 10000

    from scipy.optimize import minimize
    to_min = lambda a : np.std(energy + a*pop_f)
    res    = minimize(to_min, x0=[0.0])
    energy_corrected = energy + res.x*pop_f

    e1  = np.mean(energy)
    de1 = np.std(energy)/np.sqrt(float(len(energy)))
    e2  = np.mean(energy_corrected)
    de2 = np.std(energy_corrected)/np.sqrt(float(len(energy_corrected)))

    import matplotlib.pyplot as plt

    plt.subplot(311)
    plt.plot(pop_f)

    plt.subplot(312)
    lab = "{0} +/- {1}".format(e1, de1)
    plt.plot(energy, label=lab)
    plt.legend()

    plt.subplot(313)
    lab = "{0} +/- {1}".format(e2, de2)
    plt.plot(energy_corrected, label=lab)
    plt.legend()

    plt.show()
                
methods = {
"exp_fit"       : exp_fit,
"multi_exp_fit" : multi_exp_fit,
"average"       : average,
"reblock"       : reblocked,
"pi_weighted"   : pi_weighted,
"pop_correct"   : pop_correction
}

if len(sys.argv) < 2 or not sys.argv[1] in methods:
    err = "The first argument should be one of:"
    for m in methods:
        err += "\n    "+m
    raise Exception(err)


methods[sys.argv[1]]()
