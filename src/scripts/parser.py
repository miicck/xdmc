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
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def call_on_walkers_single(func, skip_iter, filename):
    # Call the function func(i, w, x) over walkers
    # where i is the iteration, w is the weight
    # and x is the configuration in the form
    # x = [p1, p2, p3 ... ] and p1...pN are the
    # particle position vectors

    # Loop over lines in the wavefunction file
    with open(filename) as f:
        iteration = 0
        skip = False
        for l in f:
            if l.startswith("#"):
                # Got to next iteration
                iteration += 1
                if (iteration % 1000 == 0):
                    print("    Iteration {0}".format(iteration))

                # Skip certain iterations
                skip = skip_iter(iteration)
                continue

            if skip:
                continue

            # Parse the walker
            w, x = l.split(":")
            w    = float(w)
            x    = [[float(c) for c in p.split(",")] for p in x.split(";")]

            # Call the function for this walker
            func(iteration, w, x)

def call_on_walkers(func, skip_iter=lambda i: False):
    # Call the function func(i, w, x) over walkers
    # where i is the iteration, w is the weight
    # and x is the configuration in the form
    # x = [p1, p2, p3 ... ] and p1...pN are the
    # particle position vectors

    procs = []

    # Loop over wavefunction files
    wfn_n = -1
    while True:
        wfn_n += 1
        filename = "wavefunction_{0}".format(wfn_n)
        if not os.path.isfile(filename): break
        print(filename)
        call_on_walkers_single(func, skip_iter, filename)

def bin_wavefunction():
    # This will bin the wavefunction into an
    # np array of size bins x bins x bins ...
    # with the same dimensionality as the
    # wavefunction (i.e bins^Nd).

    # Convert coords to indicies
    def indicies(c, bins, scale):
        return [int(bins*(x+scale)/(2*scale)) for x in c]

    # Will be called on every walker
    def to_call(i, w, x, bins, scale, wfn):
        c = indicies(np.array(x).flatten(), bins, scale)
        try: wfn[tuple(c)] += w
        except IndexError as e: return

    # Get the iteration range etc.
    print("Parsing wavefunction...")
    start   = int(input("Start iteration: "))
    end     = int(input("End iteration  : "))
    scale = float(input("Scale          : "))
    bins    = int(input("Bins           : "))

    # Get the dimensionallity of the wavefunction
    with open("wavefunction_0") as f:
        for line in f:
            if line.startswith("#"): continue
            particle_count = len(line.split(";"))
            spatial_dimensions = len(line.split(";")[0].split(","))
            break
    config_dimensions = particle_count * spatial_dimensions

    # Initialize the wavefunction and accumulate it
    wfn = np.zeros(shape=[bins for i in range(0, config_dimensions)])
    f = lambda i,w,x : to_call(i,w,x,bins,scale,wfn)
    call_on_walkers(f, skip_iter=lambda i: i < start or i > end)

    # Normalize the wavefunction
    norm = wfn
    while True:
        try: norm = np.trapz(norm)
        except: break
    wfn /= norm

    return { 
        "wavefunction" : wfn, 
        "scale"        : scale,
        }

def parse_evolution():
    
        # Read in our data from the evolution file
        # ignoring the first line which has y axes labels on it
        lines  = open("evolution").read().split("\n")
        y_axes = lines[0].split(",")
        data   = []
        for line in lines[1:]:
                if len(line.strip()) > 0:
                        try:
                            d = [float(i) for i in line.split(",")]
                            if np.isfinite(d).all():
                                data.append(d)
                        except Exception as e:
                            print(e)
        data = list(zip(*data))
        return [y_axes, data]

def parse_wavefunction(sys_args):

        prefix = "wavefunction"
        if "nodal_surface" in sys_args:
            sys_args.remove("nodal_surface")
            prefix = "nodal_surface"

        # Combines wavefunctions across processes
        combined = None
        wfns = []
        for f in os.listdir("."):
                if not f.startswith(prefix): continue
                wfns.append(parse_wavefunction_file(f, sys_args))

        if len(wfns) == 0:
            raise Exception("Error: no wavefunction files found!")

        wfn_sizes = [len(w) for w in wfns]
        min_size  = min(wfn_sizes)
        combined  = wfns[wfn_sizes.index(min_size)]
        for w in wfns:
            if w == combined: continue
            for i in range(0, min_size):
                combined[i].extend(w[i])

        return combined

def parse_wavefunction_file(filename, sys_args):

        # Reads a wavefunction file and returns
        # [weights, x1s, x2s ...]
        # where weights_i is the weight of the i^th walker
        # xjs_i is the position of the j^th particle in the i^th walker
        # Note that the position is read as a vector, even if the
        # system is one dimensional (i.e x -> [x] for a 1d system)
        if not os.path.isfile(filename):
                print("Error parsing wavefunction from file "+filename+" does not exist!")
                quit()

        # Read lines from wavefunction file, count iterations
        #f = open(filename)
        #lines = f.read().split("\n")
        #iter_start_lines = [i for i, l in enumerate(lines) if "#" in l]
        #total_iters = len(iter_start_lines)
        #f.close()
        total_iters = -1
        with open("evolution") as f:
            for l in f:
                total_iters += 1

        # Work out which iterations to sample
        if len(sys_args) == 1:
            
            # Read the last abs(sys_args[0]) iterations
            iter_start   = total_iters - abs(int(sys_args[0]))
            iter_end     = total_iters
            iter_spacing = 1

        elif len(sys_args) == 2:

            if int(sys_args[0]) >= 0:

                # Read start and end index from args 
                iter_start   = int(sys_args[0])
                iter_end     = int(sys_args[1])
                iter_spacing = 1

            else:

                # If we were given a negative number, sample
                # the last abs(sys_args[0]) iterations
                # (Assume second number is spacing)
                iter_start   = total_iters - abs(int(sys_args[0]))
                iter_end     = total_iters
                iter_spacing = int(sys_args[1])

        elif len(sys_args) == 3:
            
            # Read in start, end, spacing
            iter_start   = int(sys_args[0])
            iter_end     = int(sys_args[1])
            iter_spacing = int(sys_args[2])

        print("Reading {0} from iteration {1} to {2} in steps of {3}".format(
               filename, iter_start, iter_end, iter_spacing))

        data = []
        with open(filename) as f:
            iter_index = 0
            for line in f:
                if line[0] == "#":
                    iter_index += 1
                    continue
                if iter_index < iter_start: continue
                if iter_index > iter_end: break
                if iter_index % iter_spacing != 0: continue
                weight, coords = line.split(":")
                dat = [float(weight)]
                for particle in coords.split(";"):
                    try:
                        dat.append([float(w) for w in particle.split(",")])
                    except:
                        continue
                data.append(dat)

        if len(data) == 0:
            fs = "Error: iterations {0} to {1} out of range for wavefuction file {2}"
            print(fs.format(iter_start, iter_end, filename))
            quit()

        data = list(map(list, zip(*data)))
        return data

def parse_wavefunction_file_old(filename, iter_start=0, iter_end=None):
        
        # Parse wavefunction into the form
        # wavefunction[iteration number][walker number] = [weight, x0, x1...]
        # where x0, x1 ... are the coordinates of the particles
        f0 = open(filename)
        lines = f0.read().split("\n")
        f0.close()

        # Get the start location of each iteration
        i_iter_start = []
        for i, l in enumerate(lines):
                if not "#" in l: continue
                i_iter_start.append(i+1)
        i_iter_start.append(len(lines)+1)
        
        # Validate the requested number of iterations
        if iter_end is None or iter_end > len(i_iter_start) - 1:
                iter_end = len(i_iter_start)-1

        wavefunction = []
        iteration    = []

        # Read the requested number of iterations
        for i in range(i_iter_start[iter_start],
                       i_iter_start[iter_end]):

                if i >= len(lines): l = "#"
                else: l = lines[i]
                if "#" in l:
                        wavefunction.append(iteration)
                        iteration = []
                        continue

                weight, coord_data = l.split(":")
                particles = coord_data.split(";")[0:-1]

                dat = [float(weight)]
                for p in particles:
                        x = np.array([float(xi) for xi in p.split(",")])
                        dat.append(x)

                iteration.append(dat)

        return wavefunction


