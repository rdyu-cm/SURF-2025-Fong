#!/bin/bash

#SBATCH -A chm240045       # allocation name
#SBATCH --job-name=slab_dft    # Set job name
#SBATCH --partition=highmem # Set the partition 
#SBATCH --qos=cpu           # Set the QoS
#SBATCH --nodes=1            # Do not change unless you know what your doing (it set the number of nodes (do not change for non-mpi jobs))
#SBATCH --ntasks-per-node=128   # Number of MPI ranks per node (one rank per GPU)
#SBATCH --mem=900GB         # Here set to ~1TB (~22GB per core set above )[limited to 3022GB per node]
#SBATCH --time=2:00:00      # Set the max time limit 

module load gcc/11.2.0  openmpi/4.0.6 cp2k
module load conda

conda activate ase-mace12

python3 run_dft_calcs_1.py