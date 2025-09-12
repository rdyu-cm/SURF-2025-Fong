"""
Script for converting Na structures to hydronium ion systems for MACE calculations.

necessary input files:
- Na structure file locations (will need to change for specific system/locations)
    
output files:
- water and hydronium xyz files (will need to change to desired location)

notes: 
- the current implementation assumes a certain number of systems are just water
- this can be changed in the file writing section at the bottom
"""

import numpy as np
import os
import glob
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from tqdm import tqdm
import pickle
from ase.io import read, write
from ase import Atoms

target_dir = '/home/rdyu/DFT_Tutorial/dft_calcs/'

systems = []
for i in range(4):
    for j in range(5):
        systems.append(f'{target_dir}structure{j + i*50}/trajectory-input.xyz')
n_nacl = [6, 3, 3, 0, 0]
n_naoh = [0, 3, 0, 6, 0]
n_hcl  = [0, 0, 3, 0, 6]

def replace_cl_with_oh(structure, num_replacements, oh_bond_length=0.96):
    positions = structure.get_positions()
    new_atoms = []

    cl_indices = [atom.index for atom in structure if atom.symbol == 'Cl']
    if len(cl_indices) < num_replacements:
        raise ValueError(f"Not enough Cl atoms to replace. Found {len(cl_indices)}, but need {num_replacements}.")

    replacements = 0
    for idx in cl_indices:
        if replacements >= num_replacements:
            break
        cl_position = positions[idx]

        # Remove the Cl atom
        structure.pop(idx - replacements)  # adjust index for each removed atom
        replacements += 1

        # Add O atom at the Cl position
        new_atoms.append(['O', cl_position])

        # Calculate H position (along the x-axis for simplicity, adjusted by bond length)
        h_position = cl_position + np.array([oh_bond_length, 0, 0])
        new_atoms.append(['H', h_position])

    # Add the new OH groups to the structure
    for symbol, pos in new_atoms:
        structure.append(Atoms(symbol, [pos])[0])
    
    return structure

def replace_na_with_h3o(structure, num_replacements):
    # Define the relative positions of H3O atoms
    h3o_relative_positions = np.array([
        [0.0, 0.0, 0.0],  # O position
        [-0.140, 0.676, 0.895],  # H1 position
        [0.878, -0.394, -0.113],  # H2 position
        [-0.681, -0.673, -0.145]  # H3 position
    ])

    positions = structure.get_positions()
    new_atoms = []

    na_indices = [atom.index for atom in structure if atom.symbol == 'Na']
    if len(na_indices) < num_replacements:
        raise ValueError(f"Not enough Na atoms to replace. Found {len(na_indices)}, but need {num_replacements}.")

    replacements = 0
    for idx in na_indices:
        if replacements >= num_replacements:
            break
        na_position = positions[idx]

        # Remove the Na atom
        structure.pop(idx - replacements)  # adjust index for each removed atom
        replacements += 1

        # Add H3O atoms at the relative positions from the Na position
        for i, symbol in enumerate(['O', 'H', 'H', 'H']):
            new_position = na_position + h3o_relative_positions[i]
            new_atoms.append([symbol, new_position])

    # Add the new H3O groups to the structure
    for symbol, pos in new_atoms:
        structure.append(Atoms(symbol, [pos])[0])

    clean_structure = Atoms(
    numbers=structure.get_atomic_numbers(),
    positions=structure.get_positions()
    )

    clean_structure.set_cell(structure.cell)
    clean_structure.set_pbc(structure.pbc)

    return clean_structure

for ion in range(4):
    for rep in range(5):
        print(systems[ion*5 + rep])
        structure = read(systems[ion*5+rep])
        if ion == 0:
            os.makedirs(f'MACE_setup_runs/water_rep{rep}/', exist_ok=True)
            write(f'MACE_setup_runs/water_rep{rep}/init.xyz', structure)
        else:
            os.makedirs(f'MACE_setup_runs/H3O_{ion}_rep{rep}/', exist_ok=True)
            structure = replace_na_with_h3o(structure, num_replacements=ion)
            print(structure.get_positions().shape)
            print(structure.get_properties)
            print(structure.arrays.keys())
            print(structure.get_atomic_numbers().shape)
            write(f'MACE_setup_runs/H3O_{ion}_rep{rep}/init.xyz', structure)
