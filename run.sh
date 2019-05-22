CORES=$(nproc --all)
mpirun -np $CORES ./dmc  
python plot.py &
python fit_hydrogen_electrons.py &
