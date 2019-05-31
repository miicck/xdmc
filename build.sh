BASE=$(pwd)
SRC=$BASE/src
BUILD=$BASE/src/build

MPI_CC=mpicc
COMP_FLAGS="-c -O3 -g -p"
LINK_FLAGS="-p -o xdmc"

rm -r $BUILD 2> /dev/null
mkdir $BUILD
cd $BUILD

$MPI_CC $COMP_FLAGS $SRC/*.cpp

cd $BASE
$MPI_CC $LINK_FLAGS $BUILD/*.o
