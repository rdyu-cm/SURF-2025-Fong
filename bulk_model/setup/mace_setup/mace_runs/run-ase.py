'''
Copies run_1 into all MACE directories and submits a batch job to run them on HPC
'''

from ase import units
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.constraints import FixAtoms
from ase.constraints import FixedLine
from ase.md import MDLogger
from ase.io import read, write
from ase.io.extxyz import read_extxyz
from ase.io.trajectory import Trajectory
import numpy as np
import time
import os
import ssl
from ase.md.npt import NPT
# from mace.calculators.mace import MACECalculator
from mace.calculators import mace_mp

############################
# SET SIMULATION TIMES
############################
init_restart = 0.0

production_steps = 50000 # 100 ps
timestep = 1

############################
# PRODUCTION
############################


home_dir = '/anvil/scratch/x-ryu3/Bulk/setup/MACE_setup/MACE_setup_runs'
cases = ['Cs','K','Na','Li','H3O','water']
for case in cases:
    for ion in range(1,4):
        for rep in range(5):
            if case == 'water':
                if ion == 1:
                    os.system(f'cp sub_1.sl {home_dir}/{case}_rep{rep}')
                    os.system(f'cp run_1.py {home_dir}/{case}_rep{rep}')
                    os.chdir(f'{home_dir}/{case}_rep{rep}')
                    if os.path.exists('final_mace_structure.xyz'):
                        break
                    os.system(f'sbatch sub_1.sl')
                else:
                    break
            else:
                os.system(f'cp sub_1.sl {home_dir}/{case}_{ion}_rep{rep}')
                os.system(f'cp run_1.py {home_dir}/{case}_{ion}_rep{rep}')
                os.chdir(f'{home_dir}/{case}_{ion}_rep{rep}')
                if os.path.exists('final_mace_structure.xyz'):
                    break
                os.system(f'sbatch sub_1.sl')












# for case in cases:
#     for ion in range(1,4):
#         for rep in range(5):
#             # Initial configuration
#             if case == 'water':
#                 if ion == 1:
#                     os.chdir(f'/home/x-ryu3/Bulk/setup/MACE_setup/MACE_setup_runs/{case}_rep{rep}')
#                     if os.path.exists('final_mace_structure.xyz'):
#                         break
#                     init_conf = read('init.xyz')
#                 else:
#                     break
#             else:
#                 os.chdir(f'/home/x-ryu3/Bulk/setup/MACE_setup/MACE_setup_runs/{case}_{ion}_rep{rep}')
#                 if os.path.exists('final_mace_structure.xyz'):
#                     break
#                 init_conf = read('init.xyz')
            
#             init_conf.set_calculator(mace_mp(dispersion=True, device='cuda', default_dtype='float32', enable_cueq=True))
#             # init_conf.set_calculator(MACECalculator('nacl-prot_swa.model', device='cuda', default_dtype='float32',enable_cueq=True))
#             #calc = mace_mp(model='medium-mpa-0', device = 'cuda', default_dtype='float32', enable_cueq=True)
#             #init_conf.calc = calc


#             # Production run at 298K
#             MaxwellBoltzmannDistribution(init_conf, temperature_K=298)
#             dyn = Langevin(init_conf, timestep * units.fs, temperature_K=298, friction=0.0025 / units.fs)
#             def write_frame():
#                 dyn.atoms.write(f'prod-traj-{init_restart}.xyz', append=True)
#             dyn.attach(write_frame, interval=1000)
#             dyn.attach(MDLogger(dyn, init_conf, f'prod-md-{init_restart}.log', header=True, stress=False, peratom=False, mode="a"), interval=1000)
#             dyn.run(production_steps)

#             init_conf.write('final_mace_structure.xyz')

