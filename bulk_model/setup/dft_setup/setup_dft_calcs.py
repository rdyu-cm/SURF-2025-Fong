'''
This script takes trajectories and converts them into a dft input file structure

input specifications:
- traj_fn: the trajectory file in each of the specified directories
- final_structure: the structure from which the boundary conditions are taken for a given trajectory

notes:
- This is set up for a specific dft.inp file. If the dft.inp file is changed, this won't work. 
- For a more robust implementation, see the interface dft set up files
'''

import ase.io
import numpy as np
import os

# Paths containing the trajectories that we want to sample from
paths = []

home_dir = '/anvil/scratch/x-ryu3/Bulk/setup/MACE_setup/MACE_setup_runs'
cations = ['Cs','K','Na','Li','H3O']
for cat in cations:
    for ion in range(1,4):
        for rep in range(5):
            paths.append(f'{home_dir}/{cat}_{ion}_rep{rep}')
for i in range(5):
    paths.append(f'{home_dir}/water_rep{i}')

traj_fn = 'prod-traj-0.0.xyz'
final_structure = 'final_mace_structure.xyz'

structures = []

# sample 50 structures from each trajectory
n_structures = 10 
pbc = []

for i in range(len(paths)):
    traj = ase.io.read(f'{paths[i]}/{traj_fn}', ':')
    sample_indices = np.random.choice(np.arange(1,len(traj)), n_structures, replace=False)
    data_file = open(f'{paths[i]}/{final_structure}')
    
    lines = data_file.readlines()
    lattice = lines[1].split()
    
    boundary_x = lattice[0].split('\"')
    boundary_y = lattice[4]
    boundary_z = lattice[8].split('\"')
    
    for j in sample_indices:
        structures.append(traj[j])
        pbc.append([boundary_x[1],boundary_y,boundary_z[0]])

ase.io.write('structures.xyz', structures)

for j in range(len(structures)):
    os.makedirs(f'dft_calcs/structure{j}', exist_ok=True)
    ase.io.write(f'dft_calcs/structure{j}/trajectory-input.xyz', structures[j])
    
    with open('dft.inp', 'r') as dft_file:
        dft_lines = dft_file.readlines()
        dft_lines[67] = f'         A {pbc[j][0]} 0.0 0.0\n'
        dft_lines[68] = f'         B 0.0 {pbc[j][1]} 0.0\n'
        dft_lines[69] = f'         C 0.0 0.0 {pbc[j][2]}\n'
        
    with open('dft_written.inp','w') as dft_file:
        dft_file.writelines(dft_lines)
    os.system(f'cp dft_written.inp dft_calcs/structure{j}/dft.inp')
