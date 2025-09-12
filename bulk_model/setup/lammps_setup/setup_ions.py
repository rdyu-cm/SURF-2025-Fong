'''
This script uses one_run.py to set up a system of nitrate, water, and cations in lammps.
This script is only for when you have >0 ion pairs
'''

import one_run
import os
import numpy as np # type: ignore



current_dir = os.getcwd()

# run parameters
timestep = 1.0
n_npt_steps = 100000
n_nvt_steps = 1000000
temp = 298 # temperature

cation = ['Cs', 'K', 'Na', 'Li']
cation_mass_array = [132.90545, 39.0983, 22.989769, 6.94] # amu
ion_charge = 1 # magnitude of ion charge
sigma_cation = [3.884, 3.332, 2.584, 2.337] # Angstrom
eps_cation = [0.1, 0.1, 0.1, 0.16013] # kcal/mol

for cation_name, cation_mass, sigma, eps in zip(cation, cation_mass_array, sigma_cation, eps_cation):
    for ion_count in range(1,4):
        n_cation = ion_count
        n_nitrate = ion_count
        n_h2o = 64
        molecule_names = ['h2o','cation','nitrate']
        number_molecules = np.array([n_h2o, n_cation, n_nitrate],dtype=int)

        # choose initial box size
        density_guess = 1.0 # g/cm^3
        approx_box_size = ((n_h2o)/6.022e23*18.02/density_guess)**(1.0/3)/1e-8
        box_size = approx_box_size*1.25 # under-pack the box for equilibration

        # LJ interaction cutoffs
        lj_cutoff = 10.0
        coul_cutoff = lj_cutoff

        system_name = 'LAMMPS_setup_runs' # CHOOSE SYSTEM NAME

        replicates = ['0','1','2','3','4']
        for replicate in replicates:
            run_name = f'{system_name}/{cation_name}_{ion_count}_rep{replicate}/'
            
            # for each replicate generate new files
            one_run.run_setup(box_size, number_molecules, molecule_names, 
                sigma, eps, cation_mass, ion_charge,
                coul_cutoff, lj_cutoff, temp, timestep, n_npt_steps, n_nvt_steps)
            
            # copy files from the setup directory to a new directory
            os.makedirs(run_name, exist_ok=True)
            os.system('cp system.in ' + run_name)
            os.system('cp system-nvt.in ' + run_name)
            os.system('cp system.in.settings ' + run_name)
            os.system('cp system.data* ' + run_name)
            os.system('cp packed_box.xyz ' + run_name + 'init.xyz')
            os.system('cp sub-npt.sl ' + run_name)
            os.system('cp sub-nvt.sl ' + run_name)

