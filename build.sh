BASE=$(pwd)
SRC=$BASE/src
BUILD=$BASE/src/build

MPI_CC=mpicc

if [ ! -z $(which mpicc.openmpi) ]
then
	MPI_CC=mpicc.openmpi
fi

if [ ! -z $(which mpic++) ]
then
	MPI_CC=mpic++
fi

LIBS="-lstdc++ -lm"
COMP_FLAGS="-c -Wall -O3 -g -p"
LINK_FLAGS="-p -o xdmc"

rm -r $BUILD 2> /dev/null
mkdir $BUILD
cd $BUILD

echo "Building with   :" $MPI_CC $COMP_FLAGS
$MPI_CC $COMP_FLAGS $SRC/*.cpp

cd $BASE
echo "Linking with    :" $MPI_CC $LINK_FLAGS
echo "Using libraries :" $LIBS
$MPI_CC $LINK_FLAGS $BUILD/*.o $LIBS
