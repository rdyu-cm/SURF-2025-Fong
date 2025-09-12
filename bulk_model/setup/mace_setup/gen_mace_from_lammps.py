"""
Script for resizing previous LAMMPS NPT simulation data and preparing MACE setup files.

The workflow:
1. Extracts box dimensions from LAMMPS data files
2. Calculates average density and volume from NPT log files
3. Rescales structures using volume ratios
4. Converts to XYZ format with proper periodic boundary conditions
5. Organizes output files for MACE setup runs

Input files expected:
- {cation}_{ion}_rep{rep}/end_npt.data (LAMMPS data files)
- {cation}_{ion}_rep{rep}/npt.log (LAMMPS log files)
- {cation}_{ion}_rep{rep}/final_npt_structure.xyz (structure files)

Output files generated:
- MACE_setup_runs/{cation}_{ion}_rep{rep}/init.xyz (rescaled structures)
- MACE_setup_runs/{cation}_{ion}_rep{rep}/log.lammps (LAMMPS logs)

Dependencies: rescale_size module, ASE, NumPy
"""

import rescale_size
import os
from ase.io import read, write
import numpy as np

def get_box(data_file):
    with open(data_file,'r') as f:
        lines = f.readlines()
    box = []
    box.append(float(lines[11].split()[0]))
    box.append(float(lines[11].split()[1]))
    box.append(float(lines[12].split()[0]))
    box.append(float(lines[12].split()[1]))
    box.append(float(lines[13].split()[0]))
    box.append(float(lines[13].split()[1]))
    return box


def get_average_density_volume(log_file):
    with open(log_file,'r') as f:
        lines = f.readlines()
    all_density = []
    all_volume = []
    for i in range(len(lines)):
        try:
            if lines[i].split()[0] == 'Time':
                for j in range(i+91,i+102):
                    all_density.append(float(lines[j].split()[5]))
                    all_volume.append(float(lines[j].split()[4]))
        except IndexError:
            pass


    return sum(all_density)/11,sum(all_volume)/11,all_volume[-1]
    


#also fix lattice creation .lt

from_path = '/home/x-ryu3/Bulk/setup/LAMMPS_setup/LAMMPS_setup_runs/'
to_path = '/home/x-ryu3/Bulk/setup/MACE_setup/'
cations = ['Cs','K','Na','Li']
lammps = '/anvil/projects/x-chm240045/software/lammps-2Aug2023/build/lmp'


for cation in cations:
    for ion in range(1,4):
        for rep in range(5):
            box = get_box(f'{from_path}{cation}_{ion}_rep{rep}/end_npt.data')
            density, volume, final_volume = get_average_density_volume(f'{from_path}{cation}_{ion}_rep{rep}/npt.log')
            print(density,volume)

            os.makedirs(f'{to_path}MACE_setup_runs/{cation}_{ion}_rep{rep}', exist_ok=True)
            os.chdir(f'{to_path}')

            input = f'final_npt_structure.xyz'
            full_input = f'{from_path}{cation}_{ion}_rep{rep}/{input}'

            rescale_size.change_box_output_then_verify(
            lammps_exe=lammps,
            box = box,
            input_xyz=full_input,
            output_xyz="rescaled_" + input,
            x_bounds=(box[0]*(volume/final_volume)**(1/3), box[1]*(volume/final_volume)**(1/3)),
            y_bounds=(box[2]*(volume/final_volume)**(1/3), box[3]*(volume/final_volume)**(1/3)), 
            z_bounds=(box[4]*(volume/final_volume)**(1/3), box[5]*(volume/final_volume)**(1/3)),
            scale = final_volume/volume,
            ion = ion,
            cation = cation
            )
            new_box = np.array(box) * (volume/final_volume)**(1/3)
            new_box = [new_box[1] - new_box[0],new_box[3]-new_box[2],new_box[5]-new_box[4]]
            atoms = read('rescaled_final_npt_structure.xyz')
            atoms.set_cell(new_box)
            atoms.set_pbc(True)
            write('rescaled_final_npt_structure.xyz',atoms)
            os.system(f'cp rescaled_final_npt_structure.xyz MACE_setup_runs/{cation}_{ion}_rep{rep}/init.xyz')
            os.system(f'cp log.lammps MACE_setup_runs/{cation}_{ion}_rep{rep}/log.lammps')