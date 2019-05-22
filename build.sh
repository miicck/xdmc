# Build the program

# Find a working copy of mpiCC
mpiccs=("mpicpc" "mpiCC" "mpicc")
for cc in "${mpiccs[@]}"
do
	found=$(which $cc)
	if [ ! -z "$found" ]; then
		MPI_CC=$found
		break
	fi
done

# Set the mpicc flags
FLAGS_OPT="-O3"
FLAGS_DEBUG="-g -O3"

MPI_CC_FLAGS=$FLAGS_OPT
if [ "$1" == "debug" ]; then
	MPI_CC_FLAGS=$FLAGS_DEBUG
fi

echo "Building with $MPI_CC $MPI_CC_FLAGS"
$MPI_CC $MPI_CC_FLAGS -o dmc main.cpp 
