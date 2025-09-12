#!/bin/bash

#SBATCH -A chm240045-gpu       # allocation name
#SBATCH --nodes=1             # Total # of nodes 
#SBATCH --ntasks-per-node=1   # Number of MPI ranks per node (one rank per GPU)
#SBATCH --gpus-per-node=2     # Number of GPUs per node
#SBATCH --time=12:00:00        # Total run time limit (hh:mm:ss)
#SBATCH -J mace               # Job name
#SBATCH -p gpu                # Queue (partition) name

# Manage processing environment, load compilers, and applications.
module purge
module load modtree/gpu
module load cuda/11.4.2
module load cudnn/cuda-11.4_8.2
module load monitor
module list
module load conda

conda activate mace314-12

#can change hidden_irreps from between 128 and 64
#the chengucb LR implementation doesn't work with enable_cueq=True, but all other versions do work with it
#change the train and test file, as well as the location of the run_train.py file based on where you installed it
#change the model based on what model you want to use
python /anvil/scratch/x-ryu3/installations/mace/mace/cli/run_train.py \
    --name="bulk_128_lr" \
    --train_file="/anvil/scratch/x-ryu3/Bulk/train/lr_train/train.xyz" \
    --energy_key="dft_energy" \
    --forces_key="dft_forces" \
    --valid_fraction=0.05 \
    --test_file="/anvil/scratch/x-ryu3/Bulk/train/lr_train/test.xyz" \
    --config_type_weights='{"Default":1.0}' \
    --E0s="average" \
    --model="MACELES" \
    --hidden_irreps='128x0e + 128x1o' \
    --r_max=6.0 \
    --batch_size=5 \
    --max_num_epochs=250 \
    --swa \
    --start_swa=125 \
    --ema \
    --ema_decay=0.99 \
    --amsgrad \
    --restart_latest \
    --enable_cueq=True \
    --device=cuda