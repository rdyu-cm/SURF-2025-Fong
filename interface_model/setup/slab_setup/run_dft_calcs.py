'''
Runs cp2k with an already set up system
'''
 
import numpy as np
import os
import subprocess
import pandas as pd
import time 

hartree_to_ev = 27.2114
bohr_to_angstrom = 0.529177

for xy in range(1,2):
    for z in range(1,2):
        os.chdir(f'dft_calcs/structure{xy}{z}/')
        if not os.path.exists('dft.out'):
            subprocess.run(f'srun cp2k.popt dft.inp >> dft.out', shell=True, check=True)
            time.sleep(5)
        os.chdir('../../')