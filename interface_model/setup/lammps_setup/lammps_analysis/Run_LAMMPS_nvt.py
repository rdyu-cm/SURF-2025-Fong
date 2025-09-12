'''
This script runs lammps with a MLP for an NVT run sequentially (not all at once)

The current script is set up for just Cs
'''

import os
import subprocess

command = '/anvil/projects/x-chm240045/software/lammps-mace/build-ampere/lmp -k on g 1 -sf kk < system-nvt.in'

cations = ['Cs']
for cation in cations:
    for rep in range(5):
        abs_path = os.path.abspath(f'{cation}_rep{rep}')
        log_file = os.path.join(abs_path, 'output.log')
        full_cmd = f'{command} > {log_file} 2>&1'
        subprocess.run(full_cmd,cwd=abs_path, shell=True)



