BASE=$(pwd)
TEST_DIR=$BASE/src/tests
SCRIPT_DIR=$BASE/src/scripts
CORES=$(nproc --all)

for DIR in $(ls $TEST_DIR)
do
	cd $TEST_DIR/$DIR

	# Run test and profile
	mpirun -np $CORES $BASE/xdmc > /dev/null
	gprof $BASE/xdmc gmon.out > profile

	# Output basic results
	echo $DIR
	echo "    "$(cat progress_0 | grep "total time")
	echo "    "energy: $(python $SCRIPT_DIR/estimate_energy.py)
	head -n 9 profile | tail -5

	# Remove output
	mv input ..
	rm * 2> /dev/null
	mv ../input .

done
