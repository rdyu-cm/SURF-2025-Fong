'''
Loops through the set up dft folder and runs cp2k on each of them sequentially
'''

import ase.io
from ase.calculators.singlepoint import SinglePointCalculator
import ase.db
import numpy as np
import os
import subprocess
import pandas as pd
import time 

hartree_to_ev = 27.2114
bohr_to_angstrom = 0.529177

all_structures = ase.io.read('structures.xyz',':')
n_structures_total = len(all_structures)

for i in range(n_structures_total):
    print(f'structures{i}')
    os.chdir(f'dft_calcs/structure{i}/')
    if not os.path.exists('dft.out'):
        subprocess.run(f'srun cp2k.popt dft.inp >> dft.out', shell=True, check=True)
        time.sleep(5)
    os.chdir('../../')