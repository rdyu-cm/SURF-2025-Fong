#!/bin/bash

#SBATCH -A chm240045       # allocation name
#SBATCH --nodes=1             # Total # of nodes 
#SBATCH --ntasks-per-node=128   # Number of MPI ranks per node (one rank per GPU)
#SBATCH --time=12:00:00        # Total run time limit (hh:mm:ss)
#SBATCH -J npt              # Job name
#SBATCH -p wholenode                # Queue (partition) name

mpirun -np $SLURM_NTASKS /anvil/projects/x-chm240045/software/lammps-2Aug2023/build/lmp -in system.in
