'''
This script takes the test set and runs the trained model on all of the structures

input files:
- test.xyz
- your_model.model

output files:
- MACE_energies.npy: a numpy array of shape (n_structures,) containing the model's calculated energies for each structure
- MACE_forces.npy: a numpy array of shape (3*n_atoms*n_structures,) containing model's calculated flattened forces for each structure
- MACE_sizes.npy: a numpy array of shape (n_structures,) containing the model's n_atoms for each structure

notes:
- this is used in conjunction with the dft_mlp_parity.py file to create a parity plot between dft and mace values
'''

import numpy as np
from mace.calculators.mace import MACECalculator
from ase.io import read
import os

fn = '/anvil/scratch/x-ryu3/Bulk/train/lr_train/test.xyz' # test set with revPBE-D3 energies and forces

# count the number of occurrences of the word 'Lattice' in the file
with open(fn, 'r') as file:
    data = file.read()
    n_structures = data.count('Lattice')

energies = np.empty(n_structures,dtype=float)
all_forces = []
calculator = '/anvil/scratch/x-ryu3/Bulk/train/lr_train/train_new_128/bulk_128_lr_stagetwo.model'
sizes = []

for i in range(n_structures):
    print(f'Calculating structure {i}')
    conf = read(fn, i)
    conf.set_calculator(MACECalculator(calculator))#, device='cuda')) #, default_dtype='float32'
    energies[i] = conf.get_potential_energy()
    structure_forces = conf.get_forces()
    sizes.append(structure_forces.shape[0])
    all_forces.extend(structure_forces.flatten())

forces = np.array(all_forces)
sizes = np.array(sizes)
os.makedirs('mlp_ase', exist_ok= True)
np.save('mlp_ase/MACE_energies.npy', energies)
np.save('mlp_ase/MACE_forces.npy', forces)
np.save('mlp_ase/MACE_sizes.npy', sizes)


