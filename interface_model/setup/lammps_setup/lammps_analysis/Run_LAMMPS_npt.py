'''
This script runs lammps with a MLP for an NPT run sequentially (not all at once)
'''

import os
import subprocess

base_dir = '/home/x-ryu3/LAMMPS_interface/setup/LAMMPS_interface_runs'
os.chdir(base_dir)

command = '/anvil/projects/x-chm240045/software/lammps-mace/build-ampere/lmp -k on g 1 -sf kk < system.in'

cations = ['Cs', 'K', 'Na', 'Li']
for cation in cations:
    for rep in range(5):
        abs_path = os.path.abspath(f'{cation}_rep{rep}')
        log_file = os.path.join(abs_path, 'output.log')
        full_cmd = f'{command} > {log_file} 2>&1'
        subprocess.run(full_cmd,cwd=abs_path, shell=True)
        os.chdir(base_dir)



