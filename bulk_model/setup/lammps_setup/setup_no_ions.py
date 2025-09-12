'''
This script uses one_run to set up a system of pure water in lammps.
'''

import one_run
import os
import numpy as np # type: ignore

print(os.environ['PATH'])

current_dir = os.getcwd()

# run parameters
timestep = 1.0
n_npt_steps = 1000000
n_nvt_steps = 10000000
temp = 298 

# THIS IS AN EXAMPLE FOR 0 ion pairs!
# ion_count is set to one to satisfy the one_run parameters, but they aren't actually used
for ion_count in range(1):
    n_cation = ion_count
    n_nitrate = ion_count
    n_h2o = 64
    molecule_names = ['h2o']
    number_molecules = np.array([n_h2o],dtype=int)

    # choose initial box size
    density_guess = 1.0 # g/cm^3
    approx_box_size = ((n_h2o)/6.022e23*18.02/density_guess)**(1.0/3)/1e-8
    box_size = approx_box_size*1.25 # under-pack the box for equilibration

    # cation parameters (Lennard Jones parameters from Fennell et al., Dang Parameters)
    cation_mass = 22.989769
    ion_charge = 1 # magnitude of ion charge
    sigma_cation = 2.584 # Angstroms
    eps_cation = 0.1 # kcal.mol

    # LJ interaction cutoffs
    lj_cutoff = 10.0
    coul_cutoff = lj_cutoff

    system_name = 'LAMMPS_Tutorial_no_ions_' + str(ion_count) + "_" # CHOOSE SYSTEM NAME

    replicates = ['0','1','2','3','4']
    for replicate in replicates:
        run_name = f'{system_name}rep{replicate}/'
        
        # for each replicate generate new files
        one_run.run_setup(box_size, number_molecules, molecule_names, 
            sigma_cation, eps_cation, cation_mass, ion_charge,
            coul_cutoff, lj_cutoff, temp, timestep, n_npt_steps, n_nvt_steps)
        
        # copy files from the setup directory to a new directory
        os.makedirs(os.path.join(current_dir,"LAMMPS_Tutorial_runs",run_name), exist_ok=True)
        os.system('cp system.in ' + os.path.join("LAMMPS_Tutorial_runs",run_name))
        os.system('cp system-nvt.in ' + os.path.join("LAMMPS_Tutorial_runs",run_name))
        os.system('cp system.in.settings ' + os.path.join("LAMMPS_Tutorial_runs",run_name))
        os.system('cp system.data* ' + os.path.join("LAMMPS_Tutorial_runs",run_name))
        os.system('cp packed_box.xyz ' + os.path.join("LAMMPS_Tutorial_runs",run_name) + 'init.xyz')
        os.system('cp sub-npt.sl ' + os.path.join("LAMMPS_Tutorial_runs",run_name))
        os.system('cp sub-nvt.sl ' + os.path.join("LAMMPS_Tutorial_runs",run_name))

