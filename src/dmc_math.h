/*
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
*/

#ifndef __MATH__
#define __MATH__

#include <math.h>
#include <vector>

inline double fexp(double x)
{
    // The exponential used in performace
    // critical parts of the code
    return exp(x);
}

int sign(double val);
double coulomb(double q1, double q2, double r);
unsigned factorial(unsigned n);

template <class T>
class permutations
{
    // This object constructs the permutations
    // of the input vector v
    public:
        permutations(std::vector<T> a)
        {
            // Use heaps algorithm to
            // generate permutations/signs
            n     = a.size();
            nfact = factorial(n);
            row_number = 0;
            s     = new double[nfact];
            p     = new T*[nfact];
            for (unsigned i=0; i<nfact; ++i)
                p[i] = new T[n];

            double sign = 1;
            unsigned c[n];
            for(unsigned i=0; i<n; ++i)
                c[i] = 0;

            output(a, sign);

            unsigned i = 0;
            while(i < n)
            {
                if (c[i] < i)
                {
                    sign  = -sign;
                    T tmp = a[i];

                    if (i % 2 == 0)
                    {
                        // i even => swap a[0], a[i]
                        a[i] = a[0];
                        a[0] = tmp;
                    }
                    else
                    {
                        // i odd => swap a[c[i]], a[i]
                        a[i] = a[c[i]];
                        a[c[i]] = tmp;
                    }
                    output(a, sign);
                    ++c[i];
                    i = 0;
                }
                else
                {
                    c[i] = 0;
                    ++i;
                }
            }
        }

        ~permutations()
        {
            // Dealocate memory
            delete[] s;
            for (unsigned i=0; i<nfact; ++i)
                delete[] p[i];
            delete[] p;
        }

        T* operator[](unsigned i) { return p[i];  }
        unsigned size()           { return nfact; }
        unsigned elements()       { return n;     }
        double sign(unsigned i)   { return s[i];  }

    private:
        T** p;
        double*  s;
        unsigned nfact;
        unsigned n;
        unsigned row_number;

        void output(std::vector<T> a, double sign)
        {
            s[row_number] = sign;
            for (unsigned i=0; i<n; ++i)
                p[row_number][i] = a[i];
            ++ row_number;
        }
};

#endif

