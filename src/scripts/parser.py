import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def parse_evolution():
    
        # Read in our data from the evolution file
        # ignoring the first line which has y axes labels on it
        lines  = open("evolution").read().split("\n")
        y_axes = lines[0].split(",")
        data   = []
        for line in lines[1:]:
                if len(line.strip()) > 0:
                        try:
                            data.append([float(i) for i in line.split(",")])
                        except Exception as e:
                            print(e)
        data = list(zip(*data))
        return [y_axes, data]

def parse_wavefunction(iter_start=0, iter_end=None):

        # Combines wavefunctions across processes
        combined = None
        for f in os.listdir("."):
                if f.startswith("wavefunction"):
                        wf = parse_wavefunction_file(f, iter_start, iter_end)
                        if combined is None:
                            combined = wf
                        else:
                            for i, wfi in enumerate(wf):
                                combined[i].extend(wfi)
        if combined is None:
                print("Error, no wavefunction files found!")
                quit()
        return combined

def parse_wavefunction_file(filename, iter_start=0, iter_end=None):

        # Reads a wavefunction file and returns
        # [weights, x1s, x2s ...]
        # where weights_i is the weight of the i^th walker
        # xjs_i is the position of the j^th particle in the i^th walker
        # Note that the position is read as a vector, even if the
        # system is one dimensional (i.e x -> [x] for a 1d system)
        if not os.path.isfile(filename):
                print("Error parsing wavefunction from file "+filename+" does not exist!")
                quit()

        f = open(filename)
        lines = f.read().split("\n")
        f.close()
        
        if iter_end is None: 
            if iter_start < 0:
                total_iters = sum([1 for l in lines if "#" in l])
                iter_end    = total_iters - 1
                iter_start  = iter_end + iter_start
            else:
                iter_end = iter_start + 1

        iter_count = -1
        data = []
        for line in lines:
            if "#" in line:
                iter_count += 1
                continue
            if iter_count < iter_start: continue
            if iter_count >= iter_end: break

            weight, coords = line.split(":")
            dat = [float(weight)]
            for particle in coords.split(";"):
                dat.append([float(w) for w in particle.split(",")])
            data.append(dat)

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
