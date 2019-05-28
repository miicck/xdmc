
licence = """/*
    DMCPP
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
*/"""

import sys
import os
for f in sys.argv[1:]:
	
	if not os.path.isfile(f):
		continue

	# Read the lines from the source
	read  = open(f)
	lines = read.read().split("\n")
	read.close()

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
	write = open(f, "w")
	write.write(licence+"\n\n")
	for i, l in enumerate(lines):
		if i in i_licence: continue
		write.write(l+"\n")
	write.close()
