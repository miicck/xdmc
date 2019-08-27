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
                if not p["read_in"]: continue
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
            max_len = max([len(p["in_name"]) for p in params])
            fs = 'progress_file << "    " << "{0:'+str(max_len)+'} "'
            for p in params:
                f.write(ws+fs.format(p["in_name"]))
                f.write(" << params::{0}".format(p["cpp_name"]))
                f.write(r' << "\n";'+"\n")
            continue
            
        f.write(l+"\n")
