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
licence = """
    XDMC
    Copyright (C) 2019 Michael Hutcheon (email mjh261@cam.ac.uk)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    For a copy of the GNU General Public License see <https://www.gnu.org/licenses/>.
"""

import os

def set_licence_cpp_h(f):
    
    # Read the lines from the source
    with open(f) as read:
        lines = read.read().split("\n")

    # Find the Licence section if it exists
    # (including trailing whitespace)
    i_licence = []
    end_found = False
    if "/*" in lines[0]:
        for i, l in enumerate(lines):
            if end_found:
                if len(l.strip()) > 0:
                    break
            if "*/" in l:
                end_found = True
            i_licence.append(i)

    # Overwrite the file with the new licence
    with open(f, "w") as write:
        write.write("/*\n"+licence+"\n*/\n\n")
        for i, l in enumerate(lines):
            if i in i_licence: continue
            write.write(l+"\n")

    print("Set licence in : " + f)

def set_licence_python(f):
    
    # Read the lines from the source
    with open(f) as read:
        lines = read.read().split("\n")

    i_licence = []
    for i, l in enumerate(lines):
        if l.startswith("#") or len(l.strip()) == 0:
            i_licence.append(i)
        else:
            break

    with open(f, "w") as write:
        for l in licence.split("\n"):
            write.write("# "+l+"\n")
        for i, l in enumerate(lines):
            if i in i_licence: continue
            write.write(l+"\n")

    print("Set licence in : " + f)

def set_licence(f):
    if not os.path.isfile(f): return
    if f.endswith(".cpp") or f.endswith(".h"):
        set_licence_cpp_h(f)
    if f.endswith(".py"):
        set_licence_python(f)

base = "/".join(os.path.realpath(__file__).split("/")[0:-2])

for f in os.listdir(base):
    set_licence(base+"/"+f)

for f in os.listdir(base+"/gen_code"):
    set_licence(base+"/gen_code/"+f)

for f in os.listdir(base+"/scripts"):
    set_licence(base+"/scripts/"+f)

