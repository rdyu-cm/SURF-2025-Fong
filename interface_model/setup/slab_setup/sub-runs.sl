#!/bin/bash

#SBATCH -A chm240045       # allocation name
#SBATCH --nodes=1             # Total # of nodes 
#SBATCH --ntasks-per-node=128   # Number of MPI ranks per node (one rank per GPU)
#SBATCH --time=2:00:00        # Total run time limit (hh:mm:ss)
#SBATCH -J cpmd               # Job name
#SBATCH -p wholenode                # Queue (partition) name

module load gcc/11.2.0  openmpi/4.0.6 cp2k
module load conda

conda activate ase-mace12
ulimit -s unlimited

python3 run_dft_calcs_4.py