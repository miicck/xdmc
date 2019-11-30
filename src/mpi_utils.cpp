#include <mpi.h>
#include "mpi_utils.h"
#include "params.h"
#include "catch.h"

double mpi_average(double val)
{
    // Get the average of val across proccesses
    double res;
    MPI_Allreduce(&val, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    res /= double(params::np);
    return res;
}

double mpi_sum(double val)
{
    // Get the sum of val across proccesses
    double res;
    MPI_Allreduce(&val, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return res;
}

int mpi_sum(int val)
{
    // Get the sum of val across proccesses
    int res;
    MPI_Allreduce(&val, &res, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    return res;
}

// MPI unit tests
TEST_CASE("Basic MPI tests", "[mpi]")
{
    // Test MPI_SUM
    SECTION("MPI sum test")
    {
        int to_sum = 1;
        REQUIRE(mpi_sum(to_sum) == params::np);
        double to_sum_d = 1.0;
        REQUIRE(mpi_sum(to_sum_d) == double(params::np));
    }

    // Test MPI_SUM as average
    SECTION("MPI average test")
    {
        double to_average = 1.0;
        REQUIRE(mpi_average(to_average) == 1.0);
    }
}
