#!/usr/bin/python
import sys; sys.dont_write_bytecode = True # Dont generate .pyc files
import shutil
import os

# Flags controlling build
COMPILERS     = "mpic++ mpic++.openmpi mpicc mpicpc"
COMPILE_FLAGS = "-c -Wall -g -O3"
LINK_FLAGS    = "-o xdmc"
LIBS          = "-lstdc++ -lm"

# Check if clean requested
if "clean" in sys.argv:
    # Remove all build-related stuff
    # (including generated code)
    if os.path.isdir("src/build"): shutil.rmtree("src/build")
    if os.path.isfile("xdmc"): os.remove("xdmc")
    if os.path.isfile("src/params.cpp"): os.remove("src/params.cpp")
    if os.path.isfile("src/params.h"): os.remove("src/params.h")
    quit()

# Find a compiler that works
import subprocess
for c in COMPILERS.split():
    found = subprocess.check_output(["which", c])
    if len(found.strip()) > 0:
        COMPILER = c
        break

# This will generate code
import src.gen_code.gen_params

try:
    # Attempt to import multiprocessing
    from multiprocessing import cpu_count, Process
    from time import sleep

except ImportError:
    # Fallback to serial case
    cpu_count = lambda : 1

def compile_cpp(cpp):
    global COMPILER, COMPILE_FLAGS

    # Get the relative paths to the .cpp file and
    # the .o file we want to create/update
    ofile = "src/build/"+cpp.replace(".cpp",".o")
    cfile = "src/"+cpp

    # Check if we need to compile this file
    if os.path.isfile(ofile):
        if os.stat(ofile).st_mtime > os.stat(cfile).st_mtime:
            # ofile is more recent than cfile
            # => don't need to recompile
            print(ofile+" already up-to-date")
            return
        else:
            os.remove(ofile)

    # Compile a cpp file with the given filename in src/
    # to an object file in build/
    cmd = COMPILER + " " + COMPILE_FLAGS + " {0} -o {1}"
    cmd = cmd.format(cfile, ofile)
    print(cmd)
    os.system(cmd)
    if not os.path.isfile(ofile):
        raise RuntimeError("Object file {0} was not generated!".format(ofile))

print("\nCompiling on {0} cores...".format(cpu_count()))

# Create the build directory
if not os.path.exists("src/build"): os.mkdir("src/build")
    
# Get the cpp files to compile
cpp_files = [f for f in os.listdir("src/") if f.endswith(".cpp")]

# Compile the c++ files
procs = []
for cpp in cpp_files:

    # Wait for a process to become available
    while len(procs) >= cpu_count():
        for p in procs:
            if not p.is_alive():
                procs.remove(p)
        # Wait before checking procs again
        sleep(0.01) 

    if cpu_count() > 1:
        p = Process(target=compile_cpp, args=(cpp,))
        p.start()
        procs.append(p)
    else:
        compile_cpp(cpp)

# Finish the compilation
for p in procs: p.join()

# Check if .o files were created succesfully
o_files = ["src/build/"+cpp.replace(".cpp",".o") for cpp in cpp_files]
for ofile in o_files:
    if not os.path.isfile(ofile):
        raise RuntimeError("Not all object files were generated successfully!")

# Check if we need to re-link the executable
link_exe = True
if os.path.isfile("xdmc"):
    link_exe = False
    for ofile in o_files:
        if os.stat("xdmc").st_mtime < os.stat(ofile).st_mtime:
            link_exe = True
            break

print("\nLinking .o files to xdmc executable...")
if link_exe:
    # Link the object files to make the executable
    cmd = COMPILER + " " + LINK_FLAGS + " " + " ".join(o_files) + " " + LIBS
    print(cmd)
    os.system(cmd)
else:
    print("xdmc executable is up-to-date.")
