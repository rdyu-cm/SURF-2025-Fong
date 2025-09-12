'''
This script is used to loop through set up lammps systems and run them one by one
'''

import os
import subprocess

base_dir = '/home/x-ryu3/Bulk/setup/LAMMPS_setup/LAMMPS_setup_runs'
os.chdir(base_dir)

command = 'mpirun -np $SLURM_NTASKS /anvil/projects/x-chm240045/software/lammps-2Aug2023/build/lmp -in system.in'

cations = ['Cs', 'K', 'Na', 'Li']
for cation in cations:
    for ions in range(1,4):
        for rep in range(5):
            abs_path = os.path.abspath(f'{cation}_{ions}_rep{rep}')
            log_file = os.path.join(abs_path, 'output.log')
            full_cmd = f'{command} > {log_file} 2>&1'
            subprocess.run(full_cmd,cwd=abs_path, shell=True)
            os.chdir(base_dir)
