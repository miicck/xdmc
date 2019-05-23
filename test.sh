BASE=$(pwd)
TEST_DIR=$BASE/tests
CORES=$(nproc --all)

for DIR in $(ls $TEST_DIR)
do
	cd $TEST_DIR/$DIR
	echo $(pwd)
	mpirun -np $CORES $BASE/dmc  
	mv input ..
	rm *
	mv ../input .
done
