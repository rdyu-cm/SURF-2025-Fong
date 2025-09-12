'''
This script goes through EMin, NPT, and NVT with a trained MLP

input files/parameters:
- your_model.model
- NPT_steps
- prod_steps (nvt steps)
- init.xyz (your starting structure, usually generated from packmol)
- the individual EMin, NPT, and NVT parameters
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
from ase.optimize import FIRE 
import pandas as pd 
import torch 

init_restart = 0.0

prod_steps = 50000 # 50 ps
npt_length = 50000
timestep = 1
interval_npt = 1000
interval_nvt = 1000

#model location
model = '/anvil/scratch/x-ryu3/Bulk/train/lr_train/train_new_128/bulk_128_lr_stagetwo.model'

thermo_data = []

#setting the calculator, use enable_cueq=True if not using chengucb model
atoms = read('init.xyz')
calc = mace_mp(model = model, dispersion=True, device='cuda', default_dtype='float32', enable_cueq=True)
atoms.calc = calc

#EMin
print('START EMIN')
opt = FIRE(atoms, trajectory='opt.traj')
opt.run(fmax=0.2, steps=500)
atoms.write('end_min.xyz')
print('END_EMIN')

#NPT implementation, NPT dimensions will be in all directions though (mask isn't set)
print('START NPT')
MaxwellBoltzmannDistribution(atoms, temperature_K=298)
npt = NPT(atoms, timestep=units.fs, temperature_K=298, 
            externalstress=units.bar*1.01325, 
            ttime=1000*units.fs,
            pfactor=10000*units.fs,
            trajectory='npt.traj', loginterval=interval_npt)
npt.attach(MDLogger(npt, atoms, f'npt.log', header=True, stress=False, peratom=False, mode="a"), interval=interval_npt)
npt.run(npt_length)
print('END_NPT')

#Writing NPT thermo data
atoms.write('end_npt.xyz')
traj = read('npt.traj',':')
for i,frame in enumerate(traj):
    temp = frame.get_temperature() if hasattr(frame, 'get_temperature') else np.nan
    volume = frame.get_volume()
    mass = frame.get_masses().sum()
    density = mass / volume * 1.66054
    try:
        energy = frame.get_potential_energy()  # Total potential energy
    except:
        energy = np.nan  # If energy not available
    thermo_data.append({
        'step': i*interval_npt,
        'Pot_Energy' : energy,
        'temperature': temp,
        'volume': volume,
        'density': density
    })
write('all_npt.xyz', traj)
df = pd.DataFrame(thermo_data)
df.to_csv('npt_thermo.csv', index=False)


#NVT implementation
print('START NVT')
nvt = Langevin(atoms, units.fs, temperature_K=298, friction=0.0025 / units.fs)
def write_frame():
    nvt.atoms.write(f'all_nvt.xyz', append=True)
nvt.attach(write_frame, interval=1000)
traj_nvt = Trajectory('nvt.traj', 'w', atoms)
nvt.attach(traj_nvt.write, interval=1000)
nvt.attach(MDLogger(nvt, atoms, f'nvt.log', header=True, stress=False, peratom=False, mode="a"), interval=interval_nvt)
nvt.run(prod_steps)
print('END NVT')
atoms.write('end_nvt.xyz')

#writing NVT thermo data
thermo_data = []
traj = read('nvt.traj',':')
for i,frame in enumerate(traj):
    temp = frame.get_temperature() if hasattr(frame, 'get_temperature') else np.nan
    volume = frame.get_volume()
    mass = frame.get_masses().sum()
    density = mass / volume * 1.66054
    try:
        energy = frame.get_potential_energy()  # Total potential energy
    except:
        energy = np.nan  # If energy not available
    thermo_data.append({
        'step': i*interval_npt,
        'Pot_Energy' : energy,
        'temperature': temp,
        'volume': volume,
        'density': density
    })
df = pd.DataFrame(thermo_data)
df.to_csv('nvt_thermo.csv', index=False)
        

