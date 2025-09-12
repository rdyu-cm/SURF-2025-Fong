#!/bin/bash

#SBATCH -A chm240045-gpu       # allocation name
#SBATCH --nodes=1             # Total # of nodes 
#SBATCH --ntasks-per-node=1   # Number of MPI ranks per node (one rank per GPU)
#SBATCH --gpus-per-node=1     # Number of GPUs per node
#SBATCH --time=1:00:00        # Total run time limit (hh:mm:ss)
#SBATCH -J mace               # Job name
#SBATCH -p gpu                # Queue (partition) name

# Manage processing environment, load compilers, and applications.
module purge
module load modtree/gpu
#module load intel-mkl/2020.4.304 gmp/6.2.1 mpfr/3.1.6 mpc/1.1.0 zlib/1.2.11 gcc/8.4.1 numactl/2.0.14 openmpi/4.0.6 cuda/11.2.2
module load python
module load cuda/12.0.1
module list

source /anvil/projects/x-chm240045/venv/mace-venv/bin/activate

/anvil/projects/x-chm240045/software/lammps-mace/build-ampere/lmp -k on g 1 -sf kk < system.in