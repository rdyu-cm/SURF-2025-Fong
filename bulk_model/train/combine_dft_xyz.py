'''
This script creates a combined_dft.xyz file that contains all of the dft output energies
and forces, converted to the right units for training.

input parameters:
- from_dir: the directory which holds the dft structures.xyz file as well as the dft_calcs directory
'''

import ase.io
from ase.calculators.singlepoint import SinglePointCalculator
import ase.db
import numpy as np
import os
import pandas as pd
import sys

hartree_to_ev = 27.2114
bohr_to_angstrom = 0.529177

from_dir = '/anvil/scratch/x-ryu3/Bulk/setup/DFT_setup'
all_structures = ase.io.read(f'{from_dir}/structures.xyz',':')
n_structures_total = len(all_structures)

# combine all the output files (ignoring any unconverged structures)
skip = [] # from convergence check

structures = []
for i in range(n_structures_total):
    if i in skip:
        continue
    path = f'{from_dir}/dft_calcs/structure{i}'
    atoms = ase.io.read(f'{path}/trajectory-input.xyz')

    energy_file = f'{path}/ref_calc-1.ener'
    file_contents = pd.read_csv(energy_file,delim_whitespace=True,skiprows=1,
            names=['Step','Time[fs]','Kin.[a.u.]','Temp[K]','Pot.[a.u.]','Cons Qty[a.u.]','UsedTime[s]'])
    energy = float(file_contents['Pot.[a.u.]']*hartree_to_ev)

    force_file = f'{path}/forces-output.xyz'
    df = pd.read_csv(force_file,delim_whitespace=True,skiprows=2, names=['atom','x','y','z'])
    forces = np.array(df.to_numpy()[:,1:4]/bohr_to_angstrom*hartree_to_ev)

    calc = SinglePointCalculator(atoms=atoms, energy=energy, forces=forces)
    atoms.calc = calc

    atoms.set_pbc([True, True, True])

    structures.append(atoms)

ase.io.write('combined_dft.xyz', structures)