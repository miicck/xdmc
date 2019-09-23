# Set libraries, compiler flags and linker flags
# Will choose the first compiler in the list that it finds

COMPILERS="mpic++ mpicc.openmpi mpicc" 
LIBS="-lstdc++ -lm"
INCLUDE=""
COMP_FLAGS="-c -Wall -Ofast $INCLUDE"
LINK_FLAGS="-o xdmc $INCLUDE"

###############################################################
# You probably do not need to modify anything below this line #
###############################################################

# The base directory, the source directory
# and the build directory locations
BASE=$(pwd)
SRC=$BASE/src
BUILD=$BASE/src/build
GEN=$BASE/src/gen_code

# Run through compilers in preference order
# until we find one that exists on this system
for COMP in $COMPILERS
do
	if [ ! -z $(which $COMP) ]
	then
		COMPILER=$COMP
		break
	fi
done

# Check compiler has been found
if [ -z $COMPILER ]
then
    echo "None of the following compilers could be found: $COMPILERS"
    exit
fi

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

# Remove the build directory
rm -r $BUILD 2> /dev/null
