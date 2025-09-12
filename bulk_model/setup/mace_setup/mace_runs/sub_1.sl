#!/bin/bash

#SBATCH -A chm240045-gpu       # allocation name
#SBATCH --nodes=1             # Total # of nodes 
#SBATCH --ntasks-per-node=1   # Number of MPI ranks per node (one rank per GPU)
#SBATCH --gpus-per-node=1     # Number of GPUs per node
#SBATCH --time=18:00:00        # Total run time limit (hh:mm:ss)
#SBATCH -J mace               # Job name
#SBATCH -p gpu                # Queue (partition) name

# Manage processing environment, load compilers, and applications.
module purge
module load modtree/gpu
module load cuda/11.4.2
module load cudnn/cuda-11.4_8.2
module list
module load conda

conda activate ase-mace12


python3 run_1.py