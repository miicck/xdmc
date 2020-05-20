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
import os
import sys

# Add the directory with this file in to PHYTHONPATH
# and import the parameters from the params.py file
sys.path.append(os.path.dirname(os.path.realpath(__file__))) 
from params import params

# Will generate the params.h and params.cpp files from 
# params_template.h and params_template.cpp respectively

# Generate the line(s) of c++ that test the given
# named parameter according to the given test.
# These lines should evaluate to true if the parameter
# FAILS the test.
# Generates c++ with indentaton given by ws (whitespace).
def gen_test(param, test, ws=""):
    test = test.strip()

    if test == "positive":
        return ws + "if (params::{0} < 0)".format(param)

    if test.startswith("between"):

        a  = float(test.split()[1])
        b  = float(test.split()[2])
        s  = ws + "if (params::{0} < {1}\n".format(param, a)
        s += ws + " || params::{0} > {1})".format(param, b)
        return s

    if test.startswith("strings"):
        s = ws + 'if (!(params::{0} == "{1}"'.format(param, test.split()[1])
        for stest in test.split()[2:]:
            s += "\n" + ws + '   || params::{0} == "{1}"'.format(param, stest)
        s += "))"
        return s

    raise ValueError("Unkown test: "+test)

# Generate the lines of C++ to output the
# information about the allowed
# values for a given test
def gen_test_allowed(test, ws=""):
    test = test.strip()
     
    if test == "any":
        return ws + r'std::cout << "Any.\n";'+"\n";

    if test == "positive":
        return ws + r'std::cout << "Any positive value.\n";'+"\n";

    if test.startswith("between"):
        s = ws + r'std::cout << "Any value between {0} and {1} (inclusive).\n";'+"\n";
        return s.format(test.split()[1], test.split()[2])

    if test.startswith("strings"):
        s = ws      + r'std::cout << "Any of:\n";'+"\n";
        for stest in test.split()[1:]:
            s += ws + r'std::cout << "    {0}\n";'.format(stest) + "\n";
        return s

    raise ValueError("Unkown test: "+test)

# Generate the code to write the error to file if a parameter is
# not in the allowed set
def gen_test_fail_error(param, test, ws=""):
    test=test.strip()

    if test == "positive":
        return ws+r'params::error_file << "Parameter {0} must be >= 0!\n";'.format(param)

    if test.startswith("between"):
        fs = r'params::error_file << "Parameter {0} must be between {1} and {2}!\n";'
        return ws+fs.format(param, test.split()[1], test.split()[2])

    if test.startswith("strings"):
        fs = ws+r'params::error_file << "Parameter {0} must be one of the following:\n";'+"\n"
        for s in test.split()[1:]:
            fs += ws+'params::error_file << "    '+s+r'\n";'+"\n"
        return fs.format(param)

    raise ValueError("Unkown test: "+test)

# Generate the params.h file from params_template.h
def gen_params_h():

    # The warning message displayed at the top of generated files
    warning  = "// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n"
    warning += "// GENERATED FILE - DO NOT EDIT \n"
    warning += "// Generated from {0} \n"
    warning += "// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n"

    with open("params_template.h") as f: 
        lines = f.read().split("\n")

    with open("../params.h", "w") as f:

        # Write warning at top of file
        f.write(warning.format("gen_code/params_template.h"))

        for l in lines:

            # Declare parameters as extern in params.h
            if "PYTHON_GEN_PARAMS_HERE" in l:
                ws = l[0:l.find("P")]
                for p in params:
                    f.write(ws+"extern {0} {1};\n".format(p["type"], p["cpp_name"]))
                continue

            f.write(l+"\n")

# Generate the params.cpp file from params_template.cpp
def gen_params_cpp():

    # The warning message displayed at the top of generated files
    warning  = "// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n"
    warning += "// GENERATED FILE - DO NOT EDIT \n"
    warning += "// Generated from {0} \n"
    warning += "// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n"

    with open("params_template.cpp") as f:
        lines = f.read().split("\n")

    with open("../params.cpp","w") as f:

        # Write warning at op of file
        f.write(warning.format("gen_code/params_template.cpp"))

        for l in lines:

            # Define parameters in params.cpp, with default value
            if "PYTHON_GEN_PARAMS_HERE" in l:
                ws = l[0:l.find("P")]
                for p in params:
                    f.write(ws+"{0} params::{1} = {2};\n".format(p["type"],p["cpp_name"],p["default"]))
                continue

            # Parse parameters from their tag name in params.cpp
            if "PYTHON_PARSE_PARAMS_HERE" in l:
                ws = l[0:l.find("P")]
                for p in params:
                    if not "in_name" in p: continue
                    f.write(ws+"// Parse {0}\n".format(p["in_name"]))
                    f.write(ws+'if (tag == "{0}")\n'.format(p["in_name"]))
                    f.write(ws+"{\n")
                    f.write(ws+"    std::stringstream ss(split[1]);\n")
                    f.write(ws+"    ss >> params::{0};\n".format(p["cpp_name"]))
                    if not "allowed" in p:
                        f.write(ws+"    continue;\n")
                        f.write(ws+"}\n\n")
                        continue
                    f.write(gen_test(p["cpp_name"], p["allowed"], ws=ws+"    ")+"\n")
                    f.write(ws+"    {\n")
                    f.write(gen_test_fail_error(p["cpp_name"], p["allowed"], ws=ws+"        ")+"\n")
                    f.write(ws+"        return false;\n")
                    f.write(ws+"    }\n\n")
                    f.write(ws+"    continue;\n")
                    f.write(ws+"}\n\n")
                continue

            # Output params in params.cpp
            if "PYTHON_OUTPUT_PARAMS_HERE" in l:
                ws = l[0:l.find("P")]
                max_len = max([len(p["in_name"]) for p in params if "in_name" in p])
                fs = 'progress_file << "    " << "{0:'+str(max_len)+'} "'
                for p in params:
                    if not "in_name" in p: continue
                    f.write(ws+fs.format(p["in_name"]))
                    f.write(" << params::{0}".format(p["cpp_name"]))
                    f.write(r' << "\n";'+"\n")
                continue

            # Output usage info
            if "PYTHON_GEN_USAGE_INFO_HERE" in l:

                # Snippets of code to be formatted with information
                WIDTH     = 60
                DES_WIDTH = WIDTH - len("Description   : ")
                ws = l[0:l.find("P")]
                seperator    = "".join(["_" for i in range(0, WIDTH)])
                seperator    = r'std::cout << "{0}\n";'.format(seperator) + '\n' 
                title        = r'std::cout << "Parameter     = {0}\n";' + '\n'
                type_line    = r'std::cout << "Type          = {0}\n";' + '\n'
                default_line = r'std::cout << "Default value = {0}\n";' + '\n'
                description  = r'std::cout << "Description   = {0}\n";' + '\n'
                new_line     = r'std::cout << "                {0}\n";' + '\n'
                allowed_line = r'std::cout << "Allowed value(s):\n";' + '\n'
                blank_line   = r'std::cout << "\n";' + '\n'

                for p in params:
                    if not "in_name" in p: continue

                    # Write the parameter name, type and default value
                    f.write(ws+seperator)
                    f.write(ws+title.format(p["in_name"]))
                    f.write(ws+type_line.format(p["type"]))
                    f.write(ws+default_line.format(p["default"].replace('"','')))

                    # Write parameter description, keeping to a
                    # sensible line length
                    first_line = True
                    words = p["description"].split()
                    des_line = ""

                    for i, w in enumerate(words):
                        des_line += w + " "
                        if i == len(words) - 1 or len(des_line + words[i+1]) >= DES_WIDTH:
                            if first_line: 
                                f.write(ws+description.format(des_line))
                                first_line = False
                            else:
                                f.write(ws+new_line.format(des_line))
                            des_line = ""

                    f.write(ws+blank_line)
                    test = "any"
                    if "allowed" in p:
                        test = p["allowed"]
                    f.write(ws+allowed_line)
                    f.write(gen_test_allowed(test, ws=ws))

                continue
                
            f.write(l+"\n")

# Set the working directory to the gen_code dir
working_dir_before = os.getcwd()
os.chdir(os.path.dirname(os.path.realpath(__file__)))
print("\nGenerating code...")

gen_cpp = True
if os.path.isfile("../params.cpp"):
    gen_cpp = os.stat("params_template.cpp").st_mtime > os.stat("../params.cpp").st_mtime

if gen_cpp:
    print("Generating params.cpp")
    gen_params_cpp()
else:
    print("params.cpp is up-to-date")

gen_h = True
if os.path.isfile("../params.h"):
    gen_h = os.stat("params_template.h").st_mtime > os.stat("../params.h").st_mtime

if gen_h:
    print("Generating params.h")
    gen_params_h()
else:
    print("params.h is up-to-date")

# Reset the working directory
os.chdir(working_dir_before)

