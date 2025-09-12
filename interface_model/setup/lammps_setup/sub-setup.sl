#!/bin/bash

#SBATCH -A chm240045       # allocation name
#SBATCH --nodes=1             # Total # of nodes 
#SBATCH --ntasks-per-node=4   # Number of MPI ranks per node (one rank per GPU)
#SBATCH --time=1:00:00        # Total run time limit (hh:mm:ss)
#SBATCH -J compile              # Job name
#SBATCH -p standard                # Queue (partition) name

module load conda
conda activate lammps

python3 setup_piston_lammps.py
