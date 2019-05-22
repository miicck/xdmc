CORES=$(nproc --all)
mpirun -np $CORES ./dmc  
python scripts/plot_evolution.py &
python scripts/plot_wavefunction.py &
