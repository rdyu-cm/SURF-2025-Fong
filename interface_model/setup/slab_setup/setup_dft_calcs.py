'''
This script sets up a cp2k dft system for a metal slab
The current script is set up for just the slab with no previous dynamics, but the comments can be uncommented to let it work for a system with previous dynamics

input:
- home_dir: the location containing the slab.xyz
- N_KPOINTS: this can be changed by changing the xy for loop and dictates the x and y k point grid
- slab.xyz: previously created slab structure file
- boundary_y: currently set up for an orthorhombic cell, but can be changed back for a monoclinic cell (comment on the line below)

output:
- orthorhombic_slab.xyz: slab.xyz converted to be a orthorhombic cell instead of a monoclinic cell
- dft_written.inp: the dft input file with the cell and k points specified
'''
 
import ase.io
import numpy as np
import os

# Paths containing the trajectories that we want to sample from
paths = []

#For this instance, we only have one path
home_dir = '/anvil/scratch/x-ryu3/LAMMPS_interface/setup/DFT_slab_test/3pbc'
paths = [home_dir]

#traj_fn = 'prod-traj-0.0.xyz'
final_structure = 'slab.xyz'

structures = []

# sample 50 structures from each trajectory
# n_structures = 10 
pbc = []

for i in range(len(paths)):
    #traj = ase.io.read(f'{paths[i]}/{traj_fn}', ':')
    #sample_indices = np.random.choice(np.arange(1,len(traj)), n_structures, replace=False)
    with open(f'{paths[i]}/{final_structure}','r') as f:
        lines = f.readlines()
    lattice = lines[1].split()
    boundary_x = [lattice[0].split('\"')[1], lattice[1], lattice[2]]
    boundary_y = [0.0, lattice[4], lattice[5]] #change for orthorhombic cell vs monoclinic cell
    #boundary_y = [lattice[3], lattice[4], lattice[5]]
    boundary_z = [lattice[6], lattice[7], float(lattice[8].split('\"')[0])+20]
    #for j in sample_indices:
    #structures.append(traj[j])
    pbc.append([boundary_x,boundary_y,boundary_z])

#ase.io.write('structures.xyz', structures)

j = 0
for xy in range(1,2):
    for z in range(1,2):
        os.makedirs(f'dft_calcs/structure{xy}{z}', exist_ok=True)
        # ase.io.write(f'dft_calcs/structure{j}/trajectory-input.xyz', structures[j])
        atoms = ase.io.read('slab.xyz')
        atoms.set_cell([pbc[j][0][0],pbc[j][1][1],pbc[j][2][2]])
        atoms.set_pbc([True, True, True])
        ase.io.write('orthorhombic_slab.xyz',atoms)
        os.system(f'cp orthorhombic_slab.xyz dft_calcs/structure{xy}{z}/trajectory-input.xyz')        
        with open('dft.inp', 'r') as dft_file:
            with open('dft_written.inp','w') as fout:
                dft_lines = dft_file.readlines()
                for line in dft_lines:
                    new_line = line.replace('{CELL_X}',f'{pbc[j][0][0]} {pbc[j][0][1]} {pbc[j][0][2]}').replace(
                                '{CELL_Y}',f'{pbc[j][1][0]} {pbc[j][1][1]} {pbc[j][1][2]}').replace(
                                '{CELL_Z}',f'{pbc[j][2][0]} {pbc[j][2][1]} {pbc[j][2][2]}').replace(
                                '{N_KPOINTS}',str(xy)).replace(
                                '{Z_KPOINTS}',str(z))
                    fout.write(new_line)
        os.system(f'cp dft_written.inp dft_calcs/structure{xy}{z}/dft.inp')
