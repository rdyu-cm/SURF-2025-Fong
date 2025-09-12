'''
This script sets up a piston/sandwich system simulation for h3o cation

input files:
- setup_one_run_piston_lammps.py
- get_trad_lammps.py
'''

import setup_h3o_one_run_vacuum_lammps as setup_h3o_one_run_vacuum_lammps
import os
import numpy as np
import get_trad_lammps


current_dir = os.getcwd()

# run parameters
timestep = 0.001
n_steps_equil = 10000
n_steps = 50000
temp = 298 # temperature

# choose number of each species
# 1M is ~ 1 ion pair per 55 water 
n_cation = 1 #45
n_anion = 1 #45
n_h2o = 32 #2500
molecule_names = ['h2o','h3o','nitrate']
n_molecules = np.array([n_h2o, n_cation, n_anion],dtype=int)

# cation parameters
cation = ['Cs', 'K', 'Na', 'Li']
cation_mass_array = [132.90545, 39.0983, 22.989769, 6.94] # amu
ion_charge = 1 # magnitude of ion charge
sigma_cation = [3.884, 3.332, 2.584, 2.337] # Angstrom
eps_cation = [0.1, 0.1, 0.1, 0.16013] # kcal/mol

# metal parameters (these are currently set for Cu and will need to be adjusted)
slab_size = 2.0 # number of repeats of FCC unit cell
epsM = 4.72
sigmaM = 2.616

# LJ interaction cutoffs
lj_cutoff = 10.0
coul_cutoff = lj_cutoff
system_name = 'LAMMPS_interface_runs'
#for cation_name, cation_mass, sigma, eps in zip(cation, cation_mass_array, sigma_cation, eps_cation):
replicates = ['0','1','2','3','4']
sigma=0
eps=0
cation_mass=0
cation_name='h3o'
for replicate in replicates:
    run_name = f'{system_name}/h3o_rep{replicate}/'
    os.makedirs(os.path.join(current_dir,run_name), exist_ok=True)
    
    setup_h3o_one_run_vacuum_lammps.run_setup(run_name, n_molecules, molecule_names, 
        slab_size, epsM, sigmaM, 
        sigma, eps, cation_mass, cation_name, ion_charge,
        coul_cutoff, lj_cutoff, temp, timestep, n_steps, n_steps_equil)
    get_trad_lammps.remove_bonds_manual('system.data','system_nobonds.data')
    os.system('cp sub-npt.sl ' + run_name)
    os.system('cp sub-nvt.sl ' + run_name)
    os.system('cp system_nobonds.data* ' + run_name)



