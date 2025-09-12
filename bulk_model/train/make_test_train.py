'''
This script creates training and test sets from combined_dft.xyz.

input files:
- combined_dft.xyz

output files:
- test.xyz
- train.xyz

notes:
- the forces and energy labels are changed slightly to fit the mace training architecture
'''

import ase.io
import ase
import os
import random

# go from keys 'energy' and 'forces' to 'dft_energy' and 'dft_forces'
def replace_words_in_file(old_file, new_file):
    with open(old_file, 'r') as file:
        content = file.read()

    content = content.replace(':forces', ':dft_forces')
    content = content.replace(' energy', ' dft_energy')

    with open(new_file, 'w') as file:
        file.write(content)

fns = []
fns.append('/anvil/scratch/x-ryu3/Bulk/train/sr_train/combined_dft.xyz')


structures = []
for i in range(len(fns)):
    atoms = ase.io.read(fns[i], index=':')
    structures.extend(atoms)
ase.io.write('all_structures.xyz', structures)
n_test = len(structures)//10
random.shuffle(structures)
test = structures[:n_test]
train = structures[n_test:]
ase.io.write('train.xyz', train)
ase.io.write('test.xyz', test)

old_file = 'test.xyz'
replace_words_in_file(old_file, old_file)

old_file = 'train.xyz'
replace_words_in_file(old_file, old_file)