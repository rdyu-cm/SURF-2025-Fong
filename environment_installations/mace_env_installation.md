# Installing a MACE environment with CuEQ
This is for Anvil - modules may be different on other HPCs

```
conda create env_name python=3.12
conda activate env_name
conda install -c conda-forge numpy scipy matplotlib -y
pip install torch==2.6.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
conda install -c conda-forge ase -y
pip install torch-dftd
pip install cuequivariance cuequivariance-torch
pip install cuequivariance-ops-torch-cu11
```

If using the most recent version of MACE:
```
git clone https://github.com/ACEsuit/mace.git
cd mace
pip install -e .
```

If using ChengUCB version of MACE:
```
git clone https://github.com/ChengUCB/mace.git
cd mace
pip install -e .
cd ..
git clone https://github.com/ChengUCB/les.git
cd les
pip install -e .
```

Anvil Module Loading:
module load cuda/11.4.2
module load cudnn/cuda-11.4_8.2
module load conda
conda activate env_name