'''
This script is run after DFT and checks if any of the runs didn't converge
'''

import os
import numpy as np
import ase.io

def check_convergence(filename):
    text = 'SCF run converged'
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            if text in line:
                return True
    return False

all_structures = ase.io.read('structures.xyz',':')
n_structures_total = len(all_structures)

for i in range(n_structures_total):
    fn = f'dft_calcs/structure{i}/dft.out'
    converged = check_convergence(fn)
    if not converged:
        print(f'Structure {i} not converged')
