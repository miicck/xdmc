# The base directory, the source directory
# and the build directory locations
BASE=$(pwd)
SRC=$BASE/src
BUILD=$BASE/src/build

# Run through compilers in preference order
# until we find one that exists on this system
for COMP in mpic++ mpicc.openmpi mpicc
do
	if [ ! -z $(which $COMP) ]
	then
		COMPILER=$COMP
		break
	fi
done

# Set libraries, compiler flags and linker flags
LIBS="-lstdc++ -lm"
COMP_FLAGS="-c -Wall -O3 -g -p"
LINK_FLAGS="-p -o xdmc"

# Create the build directory
rm -r $BUILD 2> /dev/null
mkdir $BUILD 2> /dev/null

# Compile the .cpp files
cd $BUILD
echo "Building with   :" $COMPILER $COMP_FLAGS
$COMPILER $COMP_FLAGS $SRC/*.cpp

# Link the resulting .o files
cd $BASE
echo "Linking with    :" $COMPILER $LINK_FLAGS
echo "Using libraries :" $LIBS
$COMPILER $LINK_FLAGS $BUILD/*.o $LIBS

# Remove the build directory
rm $BUILD 2> /dev/null
