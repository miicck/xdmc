# The base directory, the source directory
# and the build directory locations
BASE=$(pwd)
SRC=$BASE/src
BUILD=$BASE/src/build
GEN=$BASE/src/gen_code

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

# Generate code
echo "Generating code ..."
cd $GEN
python gen_params.py

# Compile the .cpp files
cd $BUILD
echo "Building with   :" $COMPILER $COMP_FLAGS
$COMPILER $COMP_FLAGS $SRC/*.cpp

# Link the resulting .o files
cd $BASE
echo "Linking with    :" $COMPILER $LINK_FLAGS
echo "Using libraries :" $LIBS
$COMPILER $LINK_FLAGS $BUILD/*.o $LIBS

# Remove generated code
cd $GEN
./remove_generated.sh

# Remove the build directory
rm -r $BUILD 2> /dev/null
