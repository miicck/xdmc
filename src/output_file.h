/*
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
*/

#ifndef __OUTPUT_FILE__
#define __OUTPUT_FILE__

// Class to help with I/O. Leaves the file
// closed until it's needed.
class output_file
{
public:
    output_file() { filename = "/dev/null"; }
    output_file(std::string fn) { filename = fn; }

    template<class T>
    output_file& operator<<(T t)
    {
        if (!file.is_open()) file.open(filename, std::ofstream::trunc);
        file << t;
        if (auto_flush) flush();
        return (*this);
    }

    void open(std::string fn) { filename = fn; }
    void close() { if(file.is_open()) file.close(); }
    void flush() { if(file.is_open()) file.flush(); }

    // Set to true to automatically flush the file after each write
    bool auto_flush = false;

private:
    std::string filename;
    std::ofstream file;
};

#endif


