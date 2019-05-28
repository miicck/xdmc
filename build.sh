BASE=$(pwd)
SRC=$BASE/src
BUILD=$BASE/src/build

MPI_CC=mpicc
COMP_FLAGS="-c -O3 -g"
LINK_FLAGS="-o dmc"

rm -r $BUILD
mkdir $BUILD
cd $BUILD

$MPI_CC $COMP_FLAGS $SRC/*.cpp

cd $BASE
$MPI_CC $LINK_FLAGS $BUILD/*.o
