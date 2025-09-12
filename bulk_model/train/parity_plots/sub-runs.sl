#!/bin/bash

#SBATCH -A chm240045-gpu       # allocation name
#SBATCH --nodes=1             # Total # of nodes 
#SBATCH --ntasks-per-node=1   # Number of MPI ranks per node (one rank per GPU)
#SBATCH --gpus-per-node=1     # Number of GPUs per node
#SBATCH --time=1:00:00        # Total run time limit (hh:mm:ss)
#SBATCH -J mace               # Job name
#SBATCH -p gpu                # Queue (partition) name
#SBATCH --exclude=g008

# Manage processing environment, load compilers, and applications.
module purge
module load modtree/gpu
module load conda
module list

conda activate mace314-12

python3 run_ase_mlp.py
