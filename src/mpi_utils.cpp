#include <mpi.h>
#include "mpi_utils.h"
#include "params.h"

double mpi_average(double val)
{
    // Get the average of val across proccesses on pid 0
    double res;
    MPI_Reduce(&val, &res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    res /= double(params::np);
    return res;
}

double mpi_sum(double val)
{
    // Get the sum of val across proccesses on pid 0
    double res;
    MPI_Reduce(&val, &res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    return res;
}

int mpi_sum(int val)
{
    // Get the sum of val across proccesses on pid 0
    int res;
    MPI_Reduce(&val, &res, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return res;
}
