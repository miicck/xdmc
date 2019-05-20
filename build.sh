# Build the program

# Find a working copy of mpiCC
mpiccs=("mpiCC" "mpicc")
for cc in "${mpiccs[@]}"
do
	found=$(which $cc)
	if [ ! -z "$found" ]; then
		MPI_CC=$found
		break
	fi
done

$MPI_CC -O3 -o dmc main.cpp 
