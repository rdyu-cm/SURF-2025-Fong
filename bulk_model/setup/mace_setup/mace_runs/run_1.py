'''
Runs MACE foundation model for a single system

Usage:
- copy into each directory - done by run-ase.py
- run-ase.py will loop through all files and submit them on HPC with sub_1.sl
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

init_restart = 0.0

production_steps = 50000 # 100 ps
timestep = 1

init_conf = read('init.xyz')
init_conf.set_calculator(mace_mp(dispersion=True, device='cuda', default_dtype='float32', enable_cueq=True))
MaxwellBoltzmannDistribution(init_conf, temperature_K=298)
dyn = Langevin(init_conf, timestep * units.fs, temperature_K=298, friction=0.0025 / units.fs)
def write_frame():
    dyn.atoms.write(f'prod-traj-{init_restart}.xyz', append=True)
dyn.attach(write_frame, interval=1000)
dyn.attach(MDLogger(dyn, init_conf, f'prod-md-{init_restart}.log', header=True, stress=False, peratom=False, mode="a"), interval=1000)
dyn.run(production_steps)

init_conf.write('final_mace_structure.xyz')