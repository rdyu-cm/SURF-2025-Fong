'''
This script sets up a piston/sandwich system simulation for alkali cations in ase

input files:
- setup_one_run_piston_lammps.py
- get_trad_lammps.py
'''

import setup_one_run_piston_ase as setup_one_run_piston_ase
import os
import numpy as np

current_dir = os.getcwd()

# run parameters
timestep = 1.0
n_steps_equil = 1000
n_steps = 10000
temp = 298 # temperature

# choose number of each species
# 1M is ~ 1 ion pair per 55 water 
n_cation = 2 #45
n_anion = 2 #45
n_h2o = 110 #2500
molecule_names = ['h2o','cation','nitrate']
n_molecules = np.array([n_h2o, n_cation, n_anion],dtype=int)

# cation parameters
cation = ['Cs', 'K', 'Na', 'Li']
cation_mass_array = [132.90545, 39.0983, 22.989769, 6.94] # amu
ion_charge = 1 # magnitude of ion charge
sigma_cation = [3.884, 3.332, 2.584, 2.337] # Angstrom
eps_cation = [0.1, 0.1, 0.1, 0.16013] # kcal/mol

# metal parameters (these are currently set for Cu and will need to be adjusted)
slab_size = 3.0 # number of repeats of FCC unit cell
epsM = 4.72
sigmaM = 3.128#instead of sigma, using a the cell side length #2.616

# LJ interaction cutoffs
lj_cutoff = 10.0
coul_cutoff = lj_cutoff
system_name = 'ASE_interface_runs'
for cation_name, cation_mass, sigma, eps in zip(cation, cation_mass_array, sigma_cation, eps_cation):
    replicates = ['0','1','2','3','4']
    for replicate in replicates:
        run_name = f'{system_name}/{cation_name}_rep{replicate}/'
        os.makedirs(os.path.join(current_dir,run_name), exist_ok=True)
        setup_one_run_piston_ase.run_setup(run_name, n_molecules, molecule_names, 
            slab_size, epsM, sigmaM, 
            sigma, eps, cation_mass, cation_name, ion_charge,
            coul_cutoff, lj_cutoff, temp, timestep, n_steps, n_steps_equil)



