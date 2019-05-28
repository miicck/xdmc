BASE=$(pwd)
TEST_DIR=$BASE/src/tests
SCRIPT_DIR=$BASE/src/scripts
CORES=$(nproc --all)

for DIR in $(ls $TEST_DIR)
do
	cd $TEST_DIR/$DIR
	mpirun -np $CORES $BASE/dmc > /dev/null
	echo $DIR energy: $(python $SCRIPT_DIR/estimate_energy.py)
	mv input ..
	rm *
	mv ../input .
done
