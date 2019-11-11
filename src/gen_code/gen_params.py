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
from params import params
# Will generate the params.h and params.cpp files from 
# params_template.h and params_template.cpp respectively

# The warning message displayed at the top of generated files
warning  = "// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n"
warning += "// GENERATED FILE - DO NOT EDIT \n"
warning += "// Generated from {0} \n"
warning += "// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n"

# Generate the params.h file from params_template.h
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
            continue
            
        f.write(l+"\n")



