BASE=$(pwd)
TEST_DIR=$BASE/src/tests
SCRIPT_DIR=$BASE/src/scripts
CORES=$(nproc --all)

for DIR in $(ls $TEST_DIR)
do
	cd $TEST_DIR/$DIR

	# Run test and profile
	nice -15 mpirun -np $CORES $BASE/xdmc > /dev/null
	gprof $BASE/xdmc gmon.out > profile

	# Output basic results
	echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	echo $DIR
	echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	echo "    "$(cat progress | grep "total time")
	echo "    "energy: $(python $SCRIPT_DIR/estimate_energy.py)
	echo ""
	echo "Profiling info:"
	head -n 9 profile | tail -6
	echo ""
	echo ""

	# Remove output
	mv input ..
	rm * 2> /dev/null
	mv ../input .

done
