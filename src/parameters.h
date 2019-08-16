#ifndef __PARAMETERS__
#define __PARAMETERS__

#include <vector>
#include <fstream>
#include <sstream>
#include <map>

#include "particle.h"
#include "potential.h"
#include "simulation.h"

/*
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
        if (!file.is_open()) file.open(filename);
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
*/

class parameter_base
{
public:
    virtual std::string description() const=0;
};

// The parameter class allows generic description
// of parameters of different kinds, including default
// values, descriptions etc...
template <class T>
class parameter : public parameter_base
{
public:

    T value;
    T default_value()  const { return _default_value; }
    std::string name() const { return _name; }
    std::string description() const { return _description; }
    void operator= (const T r) { value = r; }

private:

    // By making the constructor private, parameters
    // can only be declared in a parameter_collection
    friend class parameter_collection;
    parameter(T default_val, std::string description) 
    { 
        value          = default_val;
        _default_value = default_val;
        _description   = description;
    }

    T _default_value;
    std::string _description;
};

// The collection of all parameters
// used by the program
class parameter_collection
{
public:
    void load(int argc, char** argv);

    template<class T>
    T get(std::string s)
    {
        return cast_param<T>(s)->value;
    }

    template<class T>
    void set(std::string s, T val)
    {
        cast_param<T>(s)->value = val;
    }

    // Output files
    output_file wavefunction_file;
    output_file evolution_file;
    output_file progress_file;
    output_file error_file;

private:

    template<class T>
    parameter<T>* cast_param(std::string s)
    {
        parameter<T>* pt = dynamic_cast<parameter<T>*>(params[s]);

        if (pt == nullptr)
        {
            // Cast failed
            std::stringstream err;
            err << "Could not convert parameter '" << s << "' ";
            err << "to type '" << typeid(T).name() << "'";
            throw std::runtime_error(err.str());
        }

        return pt;
    }

    // The actual parameters 
    std::map<std::string, parameter_base*> params;

    void read_input();         // Reads/parses the input file
    void output_sim_details(); // Outputs params to progress file
    int start_clock;           // The result of clock() called at load
};

extern parameter_collection params;

// Overload parameter operators to operate on parameter.value
// rather than the parameter object itself. The return type will
// be the same as the type of the object mutltiplying the parameter.

// Multiplication
template <class T, class U>
U operator* (const parameter<T>& l, U r) { return l.value * r; }
template <class T, class U>
U operator* (U l, const parameter<T>& r) { return l * r.value; }

// Division
template <class T, class U>
U operator/ (const parameter<T>& l, U r) { return l.value / r; }
template <class T, class U>
U operator/ (U l, const parameter<T>& r) { return l / r.value; }

// Adddition
template <class T, class U>
U operator+ (const parameter<T>& l, U r) { return l.value + r; }
template <class T, class U>
U operator+ (U l, const parameter<T>& r) { return l + r.value; }

// Subtraction
template <class T, class U>
U operator- (const parameter<T>& l, U r) { return l.value - r; }
template <class T, class U>
U operator- (U l, const parameter<T>& r) { return l - r.value; }

// Comparison
template <class T, class U>
U operator< (const parameter<T>& l, U r) { return l.value < r; }
template <class T, class U>
U operator< (U l, const parameter<T>& r) { return l < r.value; }

// Comparison
template <class T, class U>
U operator> (const parameter<T>& l, U r) { return l.value > r; }
template <class T, class U>
U operator> (U l, const parameter<T>& r) { return l > r.value; }

// Comparison
template <class T, class U>
U operator== (const parameter<T>& l, U r) { return l.value == r; }
template <class T, class U>
U operator== (U l, const parameter<T>& r) { return l == r.value; }

// Stream
template <class T, class U>
U operator<< (const parameter<T>& l, U r) { return l.value << r; }
template <class T, class U>
U operator<< (U l, const parameter<T>& r) { return l << r.value; }

#endif
